#!/usr/bin/perl 

use strict;
use Getopt::Long;
use List::Util qw ( min max sum );
use IPC::Run3;
no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {
    $bindir =
      dirname( abs_path(__FILE__) );    # this has to go in Begin block so happens at compile time
    $libdir = $bindir . '/../lib';
    $libdir = abs_path($libdir);        # collapses the bin/../lib to just lib
}
use lib $libdir;

use CXGN::Phylo::File;
use CXGN::Phylo::Parser;
use CXGN::Phylo::Overlap;
use CXGN::Phylo::Orthologger;
use CXGN::Phylo::Mrbayes;
use CXGN::Phylo::IdTaxonMap;

# use Devel::Cycle; # for finding circular refs - cause of memory leaks.
# find_cycle($test); # to find circular refs in $test (which might be an object, e.g.)

# find orthologs, with levels of support - either bootstrap support, or
# bayesian posterior probability from MrBayes.

# read in an alignment file. get an overlap object.
# repeat $n_bootstrap times:
# 1) get bootstrap overlap
# 2) get tree (using fasttree, quicktree (default))/
# 3) optionally reroot this tree in various ways (midpoint, mindl, minvar ...)
# 4) get orthologs for this tree
# keep track of how many times a pair is predicted orthologous
# report bootstrap support for each ortholog pair
# omitting those below threshold

my $seed_increment   = 1000;    # if seed is given on c.l. increment by this for each bootstrap.
my $do_set_error     = 0;       # 0 speeds up parsing by skipping many calls to set_error.

print "$0 ", join( " ", @ARGV ), "\n";

# defaults
my $input_file               = '';
my $n_bootstrap              = 0;
my $nongap_fraction          = 0.8;
my $kimura                   = 0;         # kimura distance correction
my $bootstrap_seed           = undef;
my $treefind_method          = 'NJ';
my $min_bs_support           = 0.0;
my $reroot_method            = 'mindl';
my $species_tree_filename    = undef;
my $query_species            = undef;
my $query_id_regex           = undef;
my $show_help                = 0;
my $ortholog_output_filename = undef;
my $additional_fasttree_options = '';
my $show_lls = 0;
my $analyze_mr_bayes = 0;       # if true, look for MrBayes output to analyze.

GetOptions(
    'inputfile|alignment=s' => \$input_file,    # contains ids to be analyzed. ids in first column
    'bootstraps|n_bootstrap=s' =>
      \$n_bootstrap,    # sequences in fasta format, superset of needed sequences
    'nongap=f'                => \$nongap_fraction,           #
    'seed=i'                  => \$bootstrap_seed,
    'treefind_method=s'       => \$treefind_method,
    'minsupport=f'            => \$min_bs_support,
    'reroot_method=s'         => \$reroot_method,
    'species_tree_filename=s' => \$species_tree_filename,
    'kimura'                  => \$kimura,
    'query_species=s'         => \$query_species,
    'query_regex=s'           => \$query_id_regex,
    'help'                    => \$show_help,
    'output_filename=s'       => \$ortholog_output_filename,
	   'ft_options=s' => \$additional_fasttree_options,
'show_lls' => \$show_lls, # show log likelihood from each bootstrap on its own line.
'mrbayes' => \$analyze_mr_bayes,
);

die help_message() . "\n" if ($show_help);

# print STDERR '@ARGV array after parsing out options: ', join( "\n", @ARGV ), "\n";
my $file_arg = shift @ARGV;    # the first argument (should be alignment filename if present).

if ($input_file) {
    if ($file_arg) {
        warn "Alignment specified both as argument ($file_arg), and as "
          . "-alignment option ($input_file). Using $input_file.\n";
    }
}
elsif ($file_arg) {
    $input_file = $file_arg;
}
if ($input_file) {
    print "Input file: $input_file \n";
}
else {
    warn "Alignment input file not specified as argument nor with -aligment option. ";
    print "Will read input alignment from stdin.\n";
}
my $base_output_filename;
my $output_tree_filename = $input_file;
$output_tree_filename =~ s/fasta/newick/;
$output_tree_filename .= '.newick' if ( !$output_tree_filename =~ /newick/ );

my $run_params_filename = '';

if ($input_file) {
    if ( -f $input_file ) {
        if ( !defined $ortholog_output_filename ) {    # if output filename not specified as option
            $base_output_filename = $input_file;
            $base_output_filename =~ s/[.]fasta$//;
            $ortholog_output_filename = $base_output_filename . '.orthologs';
        }
        else {                                         # ortholog output filename specified on cl.
            $base_output_filename = $ortholog_output_filename;
            $base_output_filename =~ s/[.]orthologs$//;
        }
        $output_tree_filename = $base_output_filename . '.newick';
        $run_params_filename  = $base_output_filename . '.run_params';

    }
    else {                                             # file doesn't exist -> exit;
        die "Usage example: bootstrap_ortholog.pl -input fam12345_alignment.fasta  "
          . "-bootstraps 100,20  -treefind ML  -reroot mindl \n"
          . "This will do 100 NJ bootstraps and 20 ML bootstraps, "
          . "with mindl rerooting (i.e. to give best consistency with species tree).";
    }
}
else {                                                 # will attempt to read alignment from stdin
    if ( defined $ortholog_output_filename ) {
        $run_params_filename = $ortholog_output_filename;
        $run_params_filename =~ s/orthologs$//;
        $run_params_filename .= 'run_params';
    }
    else {                                             # no output filename specified
        $ortholog_output_filename = 'ortholog_support.out';
        $run_params_filename      = 'run_params';
    }
}
my ($NJ_output_tree_filename, $ML_output_tree_filename) = ("NJ_$output_tree_filename", "ML_$output_tree_filename");
print "Output files: \n";
print "  Orthologs:  $ortholog_output_filename \n";
print "  Run parameters:  $run_params_filename \n";
print "  NJ tree newick:  $NJ_output_tree_filename \n";
print "  ML tree newick:  $ML_output_tree_filename \n" if($treefind_method eq 'ML');

my ( $n_NJ_bootstrap, $n_ML_bootstrap ) = ( 0, 0 );
if ( $n_bootstrap =~ /(\d+),(\d+)/ ) {
    $n_NJ_bootstrap = $1;
    $n_ML_bootstrap = $2;
}
elsif ( $n_bootstrap =~ /(\d+)/ ) {
    $n_NJ_bootstrap = $n_bootstrap;
    $n_ML_bootstrap = $n_bootstrap;
}
else {
    warn "n_bootstraps option: '$n_bootstrap' invalid. No bootstraps will be performed.\n";
    $n_bootstrap = 0;
}
print "Bootstraps: NJ: $n_NJ_bootstrap, ML: $n_ML_bootstrap \n";
#my $query_species = undef; # ($opt_q)? $opt_q: undef;
#my $query_id_regex = undef; # ($opt_Q)? $opt_Q: undef;

my $fasttree_cl = 'FastTree -wag -gamma -bionj -nosupport ';    #  -noml -nosupport';

if ( $min_bs_support > 1 ) {                                    # interpret as percentage if > 1.
    $min_bs_support *= 0.01;
    if ( $min_bs_support >= 1 ) {
        warn "Requested min bs support is > 1; Setting to zero.\n";
        $min_bs_support = 0.0;
    }
}
my $quicktree_distance_correction = ($kimura) ? 'kimura' : 'none';

#### Get the alignment:
my $align_string = '';
if ( $input_file and -f $input_file ) {
    open my $fh, "<$input_file";
    $align_string = join( "", <$fh> );
}
else {    # read alignment from stdin
    $align_string = join( "", <STDIN> );
}

if (0) {

    # fixes to $align_string:
    $align_string =~ s/IMGA[|]/IMGA_/g;    #pipes in id cause problem; replace '|' with '_'.
         #$align_string =~ s/(>[^|]+)[|][^\n]+\n/$1\n/; # delete from first pipe to end of line.
    $align_string =~ s/(>[^|]+)[|][^\n]*\n/$1\n/g;    # delete from first pipe to end of line.
}

#print STDERR $align_string, "\n";

# construct an overlap object.
my $overlap_obj = CXGN::Phylo::Overlap->new( $align_string, $nongap_fraction, $bootstrap_seed );

### Setup for orthologger.

print "Rerooting method: $reroot_method.\n";

# default species tree: 52-species tree:
my $species_newick =
'( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thelungiella_h[species=Thelungiella_halophila]:0.01, Thelungiella_s[species=Thelungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

# get the species tree if specified on c.l.
if ( defined $species_tree_filename ) {
    if ( -f $species_tree_filename ) {    # read species tree newick from file
        my $species_tree_file = CXGN::Phylo::File->new($species_tree_filename);
        $species_newick = $species_tree_file->get_tree_string();
    }
    else {
        die "species tree file: [$species_tree_filename]; no such file.\n";
    }
}

# my $species_tree = CXGN::Phylo::Parse_newick->new($species_newick, $do_set_error)->parse();
my $sparser = CXGN::Phylo::Parse_newick->new( $species_newick, $do_set_error );
my $species_tree = $sparser->parse( CXGN::Phylo::BasicTree->new() );

#find_cycle($sparser);
if ( !$species_tree ) {
    die "Species tree. Parse_newick->parse() failed to return a tree object. Newick string: "
      . $species_newick . "\n";
}

#print STDERR "# species tree: \n# $species_newick \n";

# find_cycle($species_tree);
# *********** done getting species tree ********************************************************

# ************* write run_params file **********************************************
my ( $W1, $W2 ) = ( "30", "30" );
my $run_params_string =
  sprintf( "%-$W1\s %-$W2\s \n", 'alignment source:', ($input_file) ? $input_file : 'STDIN' );
$run_params_string .= sprintf( "%-$W1\s %-$W2\s \n", 'nongap fraction:',       $nongap_fraction );
$run_params_string .= sprintf( "%-$W1\s %-$W2\s \n", 'bootstrap sample size:', $n_bootstrap );
$run_params_string .= sprintf( "%-$W1\s %-$W2\s \n", 'bootstrap seed:',        $bootstrap_seed );
$run_params_string .= sprintf( "%-$W1\s %-$W2\s \n", 'treefind method:',       $treefind_method );
$run_params_string .= sprintf( "%-$W1\s %-$W2\s \n", 'min bootstrap support:', $min_bs_support );
$run_params_string .= sprintf( "%-$W1\s %-$W2\s \n", 'rerooting method:',      $reroot_method );
$run_params_string .= sprintf( "%-$W1\s %-$W2\s \n", 'species tree file:', $species_tree_filename );
$run_params_string .=
  sprintf( "%-$W1\s %-$W2\s \n", 'NJ uses kimura correction?', $kimura ? 'yes' : 'no' );
$run_params_string .=
  sprintf( "%-$W1\s %-$W2\s \n", 'ortholog outputfile', $ortholog_output_filename );
$run_params_string .= "\nSpecies tree: \n" . $species_newick . "\n";
open my $fh_run_params, ">$run_params_filename";
print $fh_run_params $run_params_string;
close $fh_run_params;

# **********************************************************************************
my %idpair_orthocount_all      = ();    # includes results for both actual and bootstrap data.
my %idpair_actual_orthocountNJ = ();    # results for actual data only.
my %idpair_actual_orthocountML = ();    # results for actual data only.
my %idpair_bs_orthocountNJ     = ();    # results for ML analyzed bootstrap data.
my %idpair_bs_orthocountML     = ();    # results for ML analyzed bootstrap data.
my %idpair_orthocountMB        = ();    # Posterior prob. results from MrBayes.

# ********************** Analyze actual data: **************************
my $overlap_fasta_string = $overlap_obj->overlap_fasta_string('');
my $overlap_length       = $overlap_obj->get_overlap_length();

print "Overlap length: $overlap_length.\n";

open my $fh_overlap, ">", "overlap.fasta";
print $fh_overlap "$overlap_fasta_string \n";
close $fh_overlap;

# print "# fasttree command line: $fasttree_cl \n";

my $gene_tree_newick = run_quicktree( $overlap_fasta_string, $quicktree_distance_correction );

#    $gene_tree_newick =~ s/;\s*$//;    # this is gene tree newick on one line.

my $orthologger_obj = CXGN::Phylo::Orthologger->new(
    {
        'gene_tree_newick' => $gene_tree_newick,
        'species_tree'     => $species_tree,
        'reroot_method'    => $reroot_method,
        'query_species'    => $query_species,
        'query_id_regex'   => $query_id_regex
    }
);
print "Actual data: NJ analysis. Orthologger object constructed.\n";

my $rerooted_gene_tree_newick = $orthologger_obj->get_gene_tree()->generate_newick();
# my $stripped_rrgtn = $rerooted_gene_tree_newick;
# $stripped_rrgtn =~ s/\[sp.*?\]//g;
# print $stripped_rrgtn, "\n";
# my ($NJ_mllen_gene_tree_newick, $NJ_mllen_ll) = run_fasttree( $overlap_fasta_string, $fasttree_cl, $stripped_rrgtn, " -nome  -mllen ");
# print "NJ mllen ll: $NJ_mllen_ll \n";
print "Actual data; done finding rerooted NJ gene tree.\n"
      . "Gene tree written to $NJ_output_tree_filename.\n";
 open my $fh_tree, ">$NJ_output_tree_filename";
    print $fh_tree "# Actual data rerooted NJ gene tree newick:\n "
      . "$rerooted_gene_tree_newick\n";
    close $fh_tree;
my $orthologger_outstring = $orthologger_obj->ortholog_result_string();
store_orthologger_out( $orthologger_outstring,
    [ \%idpair_orthocount_all, \%idpair_actual_orthocountNJ ] );
$orthologger_obj->decircularize();

#find_cycle($orthologger_obj);

if ( $treefind_method eq 'ML' ) {
    my ( $gene_tree_newick, $ft_loglikelihood ) =
      run_fasttree( $overlap_fasta_string, $fasttree_cl , undef, $additional_fasttree_options);
print "ML loglikelihood: $ft_loglikelihood \n";
    my $rerooted_gene_tree_newick;
    my $orthologger_obj = CXGN::Phylo::Orthologger->new(
        {
            'gene_tree_newick' => $gene_tree_newick,
            'species_tree'     => $species_tree,
            'reroot_method'    => $reroot_method,
            'query_species'    => $query_species,
            'query_id_regex'   => $query_id_regex
        }
    );
    print "Actual data: ML analysis. Orthologger object constructed.\n";
    $rerooted_gene_tree_newick = $orthologger_obj->get_gene_tree()->generate_newick();

    print "Actual data; done finding rerooted ML gene tree.\n"
      . "Gene tree written to $ML_output_tree_filename.\n";
    open my $fh_tree, ">$ML_output_tree_filename";
    print $fh_tree "# Actual data rerooted ML gene tree newick:\n "
      . "$rerooted_gene_tree_newick\n";
    close $fh_tree;

    my $orthologger_outstring = $orthologger_obj->ortholog_result_string();

    store_orthologger_out( $orthologger_outstring,
        [ \%idpair_orthocount_all, \%idpair_actual_orthocountML ] );
    $orthologger_obj->decircularize();
}
print "Done with tree and ortholog finding for actual data.\n";

# **************** Done with actual data. ******************************

# **************** Now the bootstrap data. *****************************
for my $i_bs ( 1 .. $n_NJ_bootstrap ) {
    my $bs_overlap_fasta_string = $overlap_obj->bootstrap_overlap_fasta_string('');

    my $gene_tree_newick =
      run_quicktree( $bs_overlap_fasta_string, $quicktree_distance_correction );

    my $orthologger_obj = CXGN::Phylo::Orthologger->new(
        {
            'gene_tree_newick' => $gene_tree_newick,
            'species_tree'     => $species_tree,
            'reroot_method'    => $reroot_method,
            'query_species'    => $query_species,
            'query_id_regex'   => $query_id_regex
        }
    );

    my $orthologger_outstring = $orthologger_obj->ortholog_result_string();
    $orthologger_obj->decircularize();

    #find_cycle($orthologger_obj);
    store_orthologger_out( $orthologger_outstring,
        [ \%idpair_orthocount_all, \%idpair_bs_orthocountNJ ] );

    print "\r                          \r";
    print "NJ Bootstraps 1..$i_bs done.";
}    # end of NJ bootstraps
print "\n";

if ( $treefind_method eq 'ML' ) {    # also construct/analyze ML trees for the bootstrap
    for my $i_bs ( 1 .. $n_ML_bootstrap ) {
        my $bs_overlap_fasta_string = $overlap_obj->bootstrap_overlap_fasta_string('');
        my ( $gene_tree_newick, $bs_ft_ll ) =
          run_fasttree( $bs_overlap_fasta_string, $fasttree_cl, undef, $additional_fasttree_options );


#  print "show_lls: [$show_lls]\n";
        my $orthologger_obj = CXGN::Phylo::Orthologger->new(
            {
                'gene_tree_newick' => $gene_tree_newick,
                'species_tree'     => $species_tree,
                'reroot_method'    => $reroot_method,
                'query_species'    => $query_species,
                'query_id_regex'   => $query_id_regex
            }
        );

        my $orthologger_outstring = $orthologger_obj->ortholog_result_string();
        $orthologger_obj->decircularize();

        #find_cycle($orthologger_obj);
        store_orthologger_out( $orthologger_outstring,
            [ \%idpair_orthocount_all, \%idpair_bs_orthocountML ] );
	if($show_lls){
	  print "\n";
	}else{
        print "\r                                                    \r";
      }
        printf("ML Bootstraps 1..%8i done. LL: %10.4f", $i_bs, $bs_ft_ll); # , $bs_ft_ll_x);
    }
}
print "\n";    # end of loop over bootstraps

# *************** Done with bootstraps ***********************************************

my $sum_MB_count = 0;
if ($analyze_mr_bayes)
{   # do Bayesian analysis (using precomputed MrBayes output, if present), if not precomputed, skip.
        # do Bayesian analysis
    warn "Attempting to do Bayesian analysis.\n";
    my $alignment_nex_filename = $input_file . ".nex";
    $alignment_nex_filename =~ s/[.]fasta//;

    my $MB_tree_count = process_MB_out($alignment_nex_filename);

    foreach my $gene_tree_newick ( keys %$MB_tree_count ) {
        my $MB_count = $MB_tree_count->{$gene_tree_newick};
        $sum_MB_count += $MB_count;

        my $d = 1;    # distance to use
        $gene_tree_newick =~ s/(\S)([,)])/$1:$d$2/g;    # put in distances.
        $gene_tree_newick =~ s/([)])([,)])/$1:$d$2/g;

        my $orthologger_obj = CXGN::Phylo::Orthologger->new(
            {
                'gene_tree_newick' => $gene_tree_newick,
                'species_tree'     => $species_tree,
                'reroot_method'    => $reroot_method,
                'query_species'    => $query_species,
                'query_id_regex'   => $query_id_regex
            }
        );
        my $orthologger_outstring = $orthologger_obj->ortholog_result_string();

        store_orthologger_out( $orthologger_outstring,
            [ \%idpair_orthocount_all, \%idpair_orthocountMB ], $MB_count );
        $orthologger_obj->decircularize();
    }
}    # ********************** end of Bayesian analysis block *********************

# ********************* output **********************************
#print "Orthologs output filename: [$ortholog_output_filename] \n";
my $ortholog_output_string = '';

$n_bootstrap = max( $n_NJ_bootstrap, $n_ML_bootstrap );
my $iwidth         = length($n_bootstrap);
my $iformat        = '%' . $iwidth . 'i';
my $id1_old        = '';
my @sorted_idpairs = sort { $a cmp $b } keys %idpair_orthocount_all;
foreach my $idpair (@sorted_idpairs) {
    my $actual_countNJ =
      ( defined $idpair_actual_orthocountNJ{$idpair} )
      ? $idpair_actual_orthocountNJ{$idpair}
      : 0;
    my $actual_countML =
      ( defined $idpair_actual_orthocountML{$idpair} )
      ? $idpair_actual_orthocountML{$idpair}
      : 0;    # print STDERR $orthologger_outstring, "\n";
    my $bs_countNJ =
      ( defined $idpair_bs_orthocountNJ{$idpair} )
      ? $idpair_bs_orthocountNJ{$idpair}
      : 0;
    my $n_bootstrap_ML = ( $treefind_method eq 'ML' ) ? $n_ML_bootstrap : 0;
    my $n_rep_ML = ( $treefind_method eq 'ML' ) ? 1 : 0;
    my $bs_countML =
      ( defined $idpair_bs_orthocountML{$idpair} )
      ? $idpair_bs_orthocountML{$idpair}
      : 0;
    my $count_MB =
      ( defined $idpair_orthocountMB{$idpair} )
      ? $idpair_orthocountMB{$idpair}
      : 0;
    my $NJ_support = ( $n_NJ_bootstrap > 0 ) ? $bs_countNJ / $n_NJ_bootstrap : 1;
    my $ML_support = ( $n_ML_bootstrap > 0 ) ? $bs_countML / $n_ML_bootstrap : 1;
    next
      if (
            ( $actual_countNJ == 0 )
        and ( $actual_countML == 0 )
        and ($NJ_support)
        and ($ML_support)
        and ( $sum_MB_count > 0
            and ( $count_MB / $sum_MB_count ) < $min_bs_support )
      );
    my ( $id1, $id2 ) = split( " ", $idpair );

    if ( $id1 ne $id1_old ) {
        $ortholog_output_string .= sprintf(
"\n%-38s NJ($iformat)  ML($iformat) bootstrap: NJ($iformat)  ML($iformat)  MB($iformat)\n",
            $id1, 1, $n_rep_ML, $n_NJ_bootstrap, $n_bootstrap_ML, $sum_MB_count );
        $id1_old = $id1;
    }
    $ortholog_output_string .= sprintf(
        "    %-37s $iformat      $iformat                $iformat      $iformat    $iformat\n",
        $id2, $actual_countNJ, $actual_countML, $bs_countNJ, $bs_countML, $count_MB );
}

open my $fh_ortho, ">$ortholog_output_filename";
print $fh_ortho $ortholog_output_string;
close $fh_ortho;

# ************************ end of main ***************************

sub run_quicktree {
    my $overlap_fasta_string = shift;
    my $correction = shift || 'kimura';
    open my $fhtmp, ">tmp_overlap_fasta";
    print $fhtmp $overlap_fasta_string, "\n";
    close $fhtmp;

    my $overlap_stockholm_string = `sreformat stockholm tmp_overlap_fasta > tmp_overlap_stockholm`;
    my $newick_out;
    if ( $correction eq 'kimura' ) {
        $newick_out = `quicktree -kimura tmp_overlap_stockholm`;
    }
    else {
        $newick_out = `quicktree  tmp_overlap_stockholm`;
    }

    $newick_out =~ s/\s+//g;    # remove whitespace
    $newick_out =~ s/;$//;
    return $newick_out;
}

sub run_fasttree {
    my $overlap_fasta_string = shift;
    my $fasttree_cl          = shift;
    my $intree = shift;
    my $additional_options = shift;
    if($intree){
      open my $fh, ">tmp_intree";
      print $fh $intree, "\n";
      close $fh;
      $fasttree_cl .= "-intree tmp_intree ";
    }
    if($additional_options){
      $fasttree_cl .= " $additional_options ";  # e.g. -nome -mllen to just optimize branch lengths for given topo.
    }

#print stderr "fasttree cl: $fasttree_cl \n";
    my $fasttree_newick_out  = "ft_newick_default_output";
    my $fasttree_stderr_out  = "ft_stderr_default_output";
    run3( "$fasttree_cl", \$overlap_fasta_string, \$fasttree_newick_out, \$fasttree_stderr_out );

 #   print stderr "$fasttree_stderr_out \n";


# Gamma(20) LogLk = -5372.474 alpha = 1.562 rescaling lengths by 1.044   # parse ll out of ft stderr output.
    my $fasttree_loglikelihood = ( $fasttree_stderr_out =~
          / Gamma [(] \d+ [)] \s+ LogLk \s+ = \s+ ([-] \d+ [.] \d*) \s+ alpha/xm ) ? $1 : undef;
    $fasttree_newick_out =~ s/\s+//g;
    $fasttree_newick_out =~ s/;$//;
    return ( $fasttree_newick_out, $fasttree_loglikelihood );
}

sub store_orthologger_out {
    my $orthologger_outstring = shift;
    my $idpair_count_hrefs    = shift;
    my $weight                = shift || 1;
    my @orthologger_out       = split( "\n", $orthologger_outstring );
    foreach my $id_orthologs (@orthologger_out) {
        if ( $id_orthologs =~ /orthologs of\s+(\S+):\s+(.*)/ ) {
            my $id1 = $1;
            my @orthologs = split( " ", $2 );
            foreach my $id2 (@orthologs) {
                my $idpair = $id1 . "   " . $id2;
                foreach my $idpair_orthocount (@$idpair_count_hrefs) {
                    $idpair_orthocount->{$idpair} += $weight;
                }
            }    # loop over ortholog names in line of orthologger output
        }
    }    # loop over lines in orthologger output
}

# Bayesian
# we have set of gene tree newicks (for topologies visited by markov chain)
# each one has associated vector of hits in the various runs done.
# Want to analyze each topology with Orthologger, for each potential
# ortholog relationship, keep track of the number of hits (just sum over runs for now).

sub process_MB_out {
    my $alignment_nex_filename = shift;    # e.g. fam9877.nex
    my $burnin_frac = shift || 0.1;

    my $mrb_obj = CXGN::Phylo::Mrbayes->new( { 'alignment_nex_filename' => $alignment_nex_filename } );
    my $newick_count      = {};            # ref to hash of (newick(with ids): hit_count) pairs
                                           # check if MrBayes output files are present:
    my $MB_output_present = 1;
    my $run_filename = $alignment_nex_filename . ".run1.t";
    if ( !( -f $alignment_nex_filename and -f $run_filename ) ) {
      
        warn "MrBayes output files $alignment_nex_filename and $run_filename  not both available, no Bayesian analysis performed.\n";
        return $newick_count;
    }

    my $gen_param_hrefs =

      #load_params($alignment_nex_filename);
      $mrb_obj->retrieve_param_samples();    # read in from *.run?.p file
    my ( $gen_ntopo_hrefs, $newick_number_map, $number_newick_map ) =

      #  load_topologies($alignment_nex_filename);	#read in from * .run?.t file
      $mrb_obj->retrieve_topology_samples();

    # $topo_count is
    my ( $topo_count, $total_trees ) =

      #  count_topologies($gen_ntopo_hrefs);
      $mrb_obj->count_topologies($gen_ntopo_hrefs);
    my $distinct_newicks = scalar keys %$topo_count;
    print "Distinct topologies: $distinct_newicks\n";

    my @sorted_tree_numbers =
      sort { sum( @{ $topo_count->{$b} } ) <=> sum( @{ $topo_count->{$a} } ) }
      keys %$topo_count;    # tree numbers, sorted by occurrences (sum of runs).

    my $number_id_map =

      #  load_number_id_map($alignment_nex_filename);
      $mrb_obj->retrieve_number_id_map();
    print "n trees (sorted): ", scalar @sorted_tree_numbers, "\n";
    my %number_rank_map    = ();    # number is n for nth distinct topology visited by chain
    my %rank_number_map    = ();    # rank is 1 for topology with highest posterior prob. etc.
    my $index              = 1;
    my $total_trees_so_far = 0;
    my $sum_diff           = 0;
    my $L1_denom           = 0;
    foreach my $tree_number (@sorted_tree_numbers) {
        my $treecount  = sum @{ $topo_count->{$tree_number} };
        my @treecounts = @{ $topo_count->{$tree_number} };
        $L1_denom += $treecount;
        $sum_diff += max(@treecounts) - min(@treecounts);
        if ( $treecount % 2 == 1 ) {
            $sum_diff--;
            $L1_denom--;
        }

        #  print "treecounts: ", join(" ", @treecounts), "\n";
        $total_trees_so_far += $treecount;
        my $newick_with_numbers = $number_newick_map->{$tree_number};
        my $post_prob           = $treecount / $total_trees;
        $number_rank_map{$tree_number} = $index;
        $rank_number_map{$index}       = $tree_number;
        printf(
            "%4i %4i %4i %4i %4i  %10.5g %10.5g  %4i  %s\n",
            $index, @treecounts, $treecount, $total_trees_so_far, $post_prob,
            $total_trees_so_far / $total_trees,
            $tree_number, $newick_with_numbers
        );
        my $newick_with_ids =
          $mrb_obj->restore_ids_to_newick( $newick_with_numbers, $number_id_map );

        $newick_count->{$newick_with_ids} = $treecount;
        $index++;
    }
    foreach ( keys %$newick_count ) {
        print "count, newick: ", $newick_count->{$_}, " $_ \n";
    }

    printf( "topology post. distrib. L1 difference between runs: %8.5f \n",
        0.5 * $sum_diff / $L1_denom );
    print "total tree hits: $total_trees \n";
    print "distinct newicks: $distinct_newicks \n";

    my $j_run = 1;
    foreach my $gen_ntopo (@$gen_ntopo_hrefs) {
        my $string        = "Run: $j_run\n";
        my $trees_read_in = scalar keys %{$gen_ntopo};
        print "run, trees read in $j_run  $trees_read_in \n";
        my $n_burnin = int( $burnin_frac * $trees_read_in );
        my @sorted_generations = sort { $a <=> $b } keys %{$gen_ntopo};

        # print "nburnin: $n_burnin  ", scalar @sorted_generations, "\n";
        foreach my $i_gen ( @sorted_generations[ $n_burnin .. $trees_read_in - 1 ] ) {
            $string .=
                "$i_gen "
              . $number_rank_map{ $gen_ntopo->{$i_gen} } . " "
              . $gen_param_hrefs->[ $j_run - 1 ]->{$i_gen} . "\n";
        }    # loop over gens
             #  print "$string\n" if(1);
        open my $fhtp, ">run$j_run.tp";
        print $fhtp "$string\n";
        close $fhtp;
        $j_run++;
    }
    return $newick_count;
}

sub help_message {
    my $help_string =
"Usage: ortholog_support.pl alignment_xxx.fasta  -boots 100 -tree ML -seed 76543 -reroot midpoint \n"
      . "(options case-insensitive, -x or --x OK, truncated option names suffice if unambiguous, e.g. -tree ML -boots 100, etc.\n"
      . " -inputfile, -alignment ALIGNMENT_PATH     Input alignment in fasta format (no default).\n"
      . " -bootstraps, -n_bootstrap N_BOOTSTRAP     N bootstrap replicates; default: 1. (or e.g. 100,20 for 100 NJ, 20 ML replicates.)\n"
      . " -nongap FRACTION                          Alignment columns not used if have > this fraction of non-gap characters. Default: 0.8\n"
      . " -seed BOOTSTRAP_SEED                      RNG seed used in bootstrapping. Default: ? \n"
      . " -treefind_method TREE_METHOD              Choices: NJ (neighbor joining) only, ML (max likelihood) in addition of NJ. Default: NJ.\n"
      . " -minsupport MIN_BS_SUPPORT                Don't output orthologs with < this amount of bootstrap support. Default: 0.\n"
      . " -reroot_method RR_METHOD                  Choices: midpoint, minvar, mindl, none. Tree rerooting. Default: mindl.\n"
      . " -speciestree_file SPECIESTREE_PATH        Species tree newick file. Default: hard-coded tree with 55 plant species.\n"
      . " -kimura                                   Enable kimura correction of distances for NJ. Default: not enabled.\n"

#      . " -query_species QUERY_SPECIES              Output only orthologs of sequences from this species.\n"
#      . " -query_regex QUERY_REGEX                  Output only orthologs of sequences matching this (perl) regular expression.\n"
      . " -help                                     Output this help message and exit.\n"
      . " -output_filename OUTPUT_FILENAME          Output filename. Default: based on input filename, or 'ortholog_support.out' if reading from stdin.\n";
    return $help_string;
}
