#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
use List::Util qw (min max sum);


no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
  $libdir = $bindir . '/../lib'; 
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
#use lib $libdir;
use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
use lib $libdir;

use Phyml;
use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info';

use CXGN::Phylo::Overlap;
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

# The idea here is to get several NJ bootstrap trees, and use them as initial trees
# for FastTree to optimize (w.r.t. actual data, not the bootstrap resampled data).
# and then get the few with the best likelihoods, and optionally further optimize using Phyml.

my $gg_filename;  # required  cl parameter 
# by default do just fasttree, and no bootstrapping.
my $n_bs = 5;
#my ($do_ft, $n_ft_bs) = (1, 0);
my ($do_phyml, $n_phyml_bs) = (0, 0);
my $phyml_opt = 'tlr';

my $input_alignment_file = undef;
my $output_newick_file = undef;
# print "[$gg_filename] [$do_nj] [$do_ml] [$n_bs] [$ml_bs] \n";
# my $n_taxa                   = shift || 21;
my $species_tree_newick_file = undef;
my $reroot_method            = 'mindl';
my $nongap_fraction = 0.8;
my $min_overlap_length = 20;
my $max_overlap_length = 1000000000;
my $min_n_sequences = 0;
my $max_n_sequences = 1000000;
my $support = 1;		# default is to do support
my $additional_ft_options = '';



# At present:
# protein alignment only as input
# WAG
# for phyml, 4 rate categories plus invariant, alpha and p_inv estimated
GetOptions(
	   'gg_file=s'           => \$gg_filename, #
	   #	   'nj!'          => \$do_nj, # whether to do NJ tree construction for actual data
	   #  exclamation point means can use  -nj  and  -nonj
	   #	   'ft!'      => \$do_ft, # whether to do FastTree for actual data
	   'phyml!' => \$do_phyml, # whether to do Phyml tree construction for actual data
	   'phyml_opt=s' => \$phyml_opt,
	   #	   'support!'          => \$support, # -nosupport to turn off outputting of local branch support numbers.
	   'n_bs=i' => \$n_bs,	# number of NJ bootstraps to do
	   #	   'n_ft_bs=i' => \$n_ft_bs, # number of FastTree bootstraps to do
	   #	   'n_phyml_bs=i' => \$n_phyml_bs, # number of Phyml bootstraps to do
	   #    'ml_bs!' => \$ml_bs, # boolean to do ML bs or not (default: 0)
	   #	   'species_tree=s' => \$species_tree_newick_file,
	   #	   'reroot_method=s' => \$reroot_method, 
	   'nongap_fraction=f' => \$nongap_fraction,
	   'min_overlap_length=i' => \$min_overlap_length,
	   'max_overlap_length=i' => \$max_overlap_length,
	   'max_n_sequences=i' => \$max_n_sequences,
	   'input_alignment_file=s' => \$input_alignment_file,
	   'output_newick_file=s' => \$output_newick_file,
	   #	   'ft_options=s' => \$additional_ft_options,
	  );
if (!defined $gg_filename  or  ! -f $gg_filename) {
  # die "No gene-genome association file specified. Exiting. \n";
}
my $gg_hashref = store_gg_info($gg_filename);

my $fh_in = \*STDIN;
if (defined $input_alignment_file) {
  open $fh_in, "<", "$input_alignment_file" or die "File [$input_alignment_file] could not be opened for reading.\n";
}
my $fh_out;
if (defined $output_newick_file) {
  open $fh_out, ">", "$output_newick_file" or die "File [$output_newick_file] could not be opened for writing.\n";
}else{
 $fh_out = \*STDOUT;
}
my $species_tree = get_species_tree($species_tree_newick_file);


print STDERR # "gg_file: ", (defined $gg_filename)? $gg_filename: 'undefined', " \n", 
  "n_bs: $n_bs \n", 
  "do_phyml: $do_phyml \n",
  "min nongap_fraction: $nongap_fraction \n",
  "min_overlap_length: $min_overlap_length \n";
my $state              = 'idline'; # other possible values: fasta
my ( $qid, $fam_size, $taxa, $idline, $fasta);
my $support_string = ($support)? '' : '-nosupport'; # controls whether FastTree calculates local branch support numbers.

my ($n_ft, $n_phyml) = (0, 0);

while (<$fh_in>) {
  if ( /^Id\s+(\S+)\s+family/ ) {
    $qid = $1;
    my @cols = split( " ", $_ );
    ( $qid, $fam_size, $taxa) = @cols[ 1, 4, 5];
    $idline = $_;
    $state  = 'fasta';
    $fasta  = '';
  } elsif ( $state eq 'fasta' ) {

    if (/^\s*$/) { # blank line after sequence -> process the sequence.
      chomp $idline;
      my $string_to_print = "$idline   ";
      my $rng_seed = srand();
      if ($fasta ne '') {
	my $overlap_obj =
	  CXGN::Phylo::Overlap->new( $fasta, $nongap_fraction,
				     $rng_seed );
	my $alignment_overlap =
	  $overlap_obj->overlap_fasta_string('');
	my $overlap_length = $overlap_obj->get_overlap_length();
	#	print "id: $qid  overlap_length: $overlap_length  n_sequences: ", $overlap_obj->get_n_sequences(), "\n";
	$string_to_print .= "$overlap_length ";
	if ( $overlap_length >= $min_overlap_length and $overlap_length <= $max_overlap_length and $overlap_obj->get_n_sequences() < $max_n_sequences) { # sufficient overlap
	  my %description_lnL = ();

	  ################## run FastTree  ######################################## 
	  # standard way:
	  my $fasttree_command = "FastTree -wag -gamma -bionj -nosupport $additional_ft_options ";
	  my ($ft_newick, $ft_lnL, $ft_cputime) = run_fasttree($alignment_overlap, $fasttree_command);
	  my $type = "FT";
	  my $description = "$type  $ft_cputime  $ft_newick";
	$ft_lnL = -1e100 if($ft_lnL eq '---');
	  $description_lnL{$description} = $ft_lnL;
	  print STDERR "$type $ft_lnL;  ";


	#   # using NJ tree as initial tree (rather than bionj done within FT)
	#   my ($njft_newick, $njft_lnL, $njft_cputime) = nj_to_ft($alignment_overlap, $alignment_overlap);
	#   $type = 'NJ->FT';
	#   $description = "$type  $njft_cputime  $njft_newick";
	# $njft_lnL = -1e100 if($njft_lnL eq '---');
	#   $description_lnL{$description} = $njft_lnL;
	#   print STDERR "$type $njft_lnL;  ";


	  # using NJ bootstrap trees as initial trees.
	  for my $i_bs (1..$n_bs) {
	    my $bs_alignment_overlap =
	      $overlap_obj->bootstrap_overlap_fasta_string('');
	    my $type = 'BS' . $i_bs . 'NJ->FT ';
	    my ($bsnjft_newick, $bsnjft_lnL, $bsnjft_cputime) = nj_to_ft($bs_alignment_overlap, $alignment_overlap);
	$bsnjft_lnL = -1e100 if($bsnjft_lnL eq '---');	    
$description_lnL{"$type  $bsnjft_cputime  $bsnjft_newick"} = $bsnjft_lnL;
	    print STDERR "$type $bsnjft_lnL;  ";
	  }			# end of bootstraps loop
	  print STDERR "\n";


	  # sort by FT likelihood 
	  my @skeys = sort {$description_lnL{$b} <=> $description_lnL{$a} } keys %description_lnL; #Sort by FastTree likelihoods (best first)
#print "A: \n  ", join("\n  ", @skeys), "\n";
	  print "Id $qid  \n";
	  my $best_ft_lnL = $description_lnL{$skeys[0]};
	  my $i = 1;
	  for my $the_key (@skeys) {

	    my ($descript, $ft_cputime, $newick) = split(" ",$the_key);
	    my $ft_lnL = $description_lnL{$the_key};
#print "QQQQQQQ: $the_key \n";
	    # phyml, if requested.
	    my ($phymlobj, $phymlnewick, $phyml_lnL, $phyml_cput) = ($do_phyml and ($i <= 2))? 
	      run_phyml($alignment_overlap, $newick, 'r') : 
		(undef, '()', 0, 0);
#print "RRRRRRRR: $the_key \n";
	    my $delta_lnL = $best_ft_lnL - $ft_lnL;
	    printf("%3i %28s %10s  %11.2f  %5.3f  %8.6f  %5.2f  %12.2f  %5.1f \n",
		   $i, $qid, $descript, $ft_lnL, $delta_lnL, $delta_lnL/abs($best_ft_lnL), $ft_cputime, $phyml_lnL, $phyml_cput);
#	    print "ZZZZ: \n", $skeys[0], "\n";
	    if ($i == 1) {
               if(!defined $phymlobj){
	      print $fh_out ("Id $qid  $descript  $ft_lnL  $ft_cputime   -  0  0  \n");
   my $taxonified_newick =
	      taxonify_newick( $newick, $gg_hashref );
	my    $rr_tax_newick =
	      newick2rerootednewick( $taxonified_newick,
				     $reroot_method, $species_tree );


	      print $fh_out "$descript  $rr_tax_newick \n\n";
           }else{
	      print $fh_out ("Id $qid  $descript  $ft_lnL  $ft_cputime   r  $phyml_lnL $phyml_cput \n");
	      $phymlnewick =~ s/\s//g; # remove whitespace
	      $phymlnewick =~ s/;\s*$//; # remove final ;
	      print $fh_out "$descript  $phymlnewick \n\n";
           }
	    }
	    $i++;
	  }			# end of loop over sorted trees
#print "B: ", join("\n", @skeys), "\n";
	  if ($do_phyml) {
	  #  print "XXX: ", $skeys[0], "\n";
	    my ($descript, $ft_cputime, $newick_in) = split(" ", $skeys[0]);
	 #   print "$descript  $ft_cputime  $newick_in \n";
	    my $ft_lnL = $description_lnL{$skeys[0]};
	    my ($phymlobj, $phymlnewick, $phyml_lnL, $phyml_cput) = run_phyml($alignment_overlap, $newick_in, $phyml_opt);
	    print "Id $qid  $descript  $ft_lnL  $ft_cputime  $phyml_opt  $phyml_lnL $phyml_cput \n";
	    print $fh_out ("Id $qid  $descript  $ft_lnL  $ft_cputime  $phyml_opt  $phyml_lnL $phyml_cput \n");
	    $phymlnewick =~ s/\s//g; # remove whitespace
	    $phymlnewick =~ s/;\s*$//; # remove final ;
	    print $fh_out "$descript  $phymlnewick \n\n";
	  }
	}
      }
      $state = 'idline';
    } else {			# not an empty line - append to $fasta
      $fasta .= $_;
    }
  }
}

sub nj_to_ft{  # get a tree using NJ and use as init tree for FastTree
  # usually want bootstrap alignment as input to NJ, then real alignment as input to FastTree
  my $nj_input_alignment = shift; #
  my $ft_input_alignment = shift;
  my $nj_newick;
  my ($bs_lnL, $bs_cputime);
  if (1) {			# use NJ for the bs trees ...
    $nj_newick = run_quicktree($nj_input_alignment);
  } else { # can use FastTree (just min. evolution, no ML) for bs trees ...
    ($nj_newick, $bs_lnL, $bs_cputime) = run_fasttree($nj_input_alignment,"FastTree -wag -gamma -bionj -noml" );
  }
  my $fasttree_command = "FastTree -wag -gamma -bionj ";
  my ($ft_newick, $ft_lnL, $ft_cpu_time) = 
    run_fasttree($ft_input_alignment, $fasttree_command, $nj_newick);
  return ($ft_newick, $ft_lnL, $ft_cpu_time);
}

sub taxonify_newick {
  my $newick        = shift;
  my $seqid_species = shift;
  # done storing gg_file info in hash %seqid_species
  $newick =~ s/\s+$//;		# remove final whitespace
  my $new_newick = $newick;
  $new_newick =~ s/ ([\(,]) \s* ([^,\):]+) \s* ([,\):])/$1$2\n$3/xg
    ;				# add newline after each leaf
  my @leaf_lines = split( "\n", $new_newick );
  my $last_bit = pop @leaf_lines;
  for (@leaf_lines) {
    if (/\[species=.*\]/) { # species info already there for this leaf - leave it alone
    } else {
      / [\(,] \s* ([^\(,\):]+) \s* $/x;
      my $seq_id = $1;
      if ( exists $seqid_species->{$seq_id} ) {
	my $species = $seqid_species->{$seq_id};
	$seq_id .= '[species=' . $species . ']';
      } else {
	warn "sequence $seq_id; no corresponding species found.\n";
      }
      s/ ([\(,]) \s* ([^\(,\):]+) \s* $/$1$seq_id/x;
    }
  }
  $new_newick = join( '', @leaf_lines ) . $last_bit;
  return $new_newick;
}

sub reroot {
  my $tree             = shift;
  my $reroot_method    = shift || 'mindl';
  my $species_tree     = shift;
  my $species_tree_arg = shift; # file name of species tree newick file. overrides default species tree

  #  print "in reroot. reroot method $reroot_method \n";
  #print $tree->generate_newick(), "\n";
  my ( $new_root, $dist_above ) = ( undef, undef );
  if ( $reroot_method eq 'none' ) {
    return $tree;
  } elsif ( $reroot_method eq 'midpoint' ) {
    ( $new_root, $dist_above ) = $tree->find_midpoint();
  } elsif ( $reroot_method eq 'minvar' ) {
    ( $new_root, $dist_above ) = $tree->min_leaf_dist_variance_point();
  } elsif ( $reroot_method eq 'mindl' ) {

  
    my $spec_bit_hash = $tree->get_species_bithash($species_tree);

    $tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);

    $species_tree->get_root()
      ->recursive_set_implicit_species_bits($spec_bit_hash);

    ( $new_root, $dist_above ) = $tree->find_mindl_node($species_tree);

  } else {
    warn "reroot option $reroot_method unknown. No rerooting performed.\n";
  }
  if ( defined $new_root ) {
    $tree->reset_root_to_point_on_branch( $new_root, $dist_above );
  } else {
    warn "In reroot. \$new_root undefined, can't reroot. Tree unchanged.\n";
  }
  return $tree;
}

sub newick2rerootednewick
  { # take a taxonified newick expression, make tree, reroot is by mindl wrt $species_tree

    # return newick expression
    my $taxonified_newick = shift;
    $taxonified_newick = put_min_branchlengths_into_newick($taxonified_newick, 0.01);
    my $reroot_method     = shift;
    my $species_tree      = shift;
    my $parser = CXGN::Phylo::Parse_newick->new( $taxonified_newick, 0 );
    my $tree   = $parser->parse( CXGN::Phylo::BasicTree->new() );
    if (!defined $tree) {
      return 'warning: Couldnt parse taxonified tree in newick2rerootnewick.';
    }
    $tree->impose_branch_length_minimum();
    $tree->show_newick_attribute("species");
    $tree->set_show_standard_species(0);
    $tree->get_root()->recursive_implicit_names();
    $tree->get_root()->recursive_implicit_species();

    # if species is empty, (because no [species=x] in newick), get it from the name...
    #$tree->set_missing_species_from_names();

    $tree->set_show_standard_species(1);
    $tree->set_species_standardizer( CXGN::Phylo::Species_name_map->new() );

    ### make the tree binary
    my @new_root_point;
    { # if root has > 2 children, reroot to pt on a root-child branch, so new root has 2 children.
      my @root_children = $tree->get_root()->get_children();
      if ( scalar @root_children > 2 ) {
	@new_root_point = (
			   $root_children[0], 0.9 * $root_children[0]->get_branch_length()
			  );
	$tree->reset_root_to_point_on_branch(@new_root_point);
      }
    }

    # binarify every non-binary node. At present doesn't attempt to choose in a smart way
    # among the various possible resolutions
    $tree->make_binary( $tree->get_root() ); # urec requires binary tree.
  
    $tree = reroot( $tree, $reroot_method, $species_tree );

    my $rerooted_tree_newick = $tree->generate_newick();

    $tree->decircularize();    # so can be properly garbage collected!
    return $rerooted_tree_newick;
  }

sub get_species_tree {
  my $species_tree_arg = shift
    ; # file name of species tree newick file. default species tree as specified below
  my $species_tree_newick;
  if ( defined $species_tree_arg ) {
    my $species_tree_file = CXGN::Phylo::File->new($species_tree_arg);
    find_cycle($species_tree_file);
    $species_tree_newick = $species_tree_file->get_tree_string();
  } else {
    $species_tree_newick = # includes most of 50 species in Sept. 2014 analysis. (but no phaseolus, phyllostachys, fraxinus, sesamum (
      '( ((ostreococcus_t[species=Ostreococcus_tauri]:1, ostreococcus_l[species=Ostreococcus_lucimarinus]:1):1, (chlamydomonas[species=Chlamydomonas_reinhardtii]:1, volvox[species=Volvox_carteri]:1):1):1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( (norway_spruce[species=Picea_abies]:1, loblolly_pine[species=Pinus_taeda]:1):1, ( amborella[species=Amborella_trichopoda]:1, ( (duckweed[species=Spirodela_polyrhiza]:1, ( date_palm[species=Phoenix_dactylifera]:1, (banana[species=Musa_acuminata]:1,( ( (switchgrass[species=Panicum_virgatum]:1, foxtail_millet[species=Setaria_italica]:1):1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1 ):1):1, ( columbine[species=Aquilegia_coerulea]:1, ( nelumbo[species=Nelumbo_nucifera]:1, ( ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, (monkey_flower[species=Mimulus_guttatus]:1, bladderwort[species=Utricularia_gibba]:1):1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, (carnation[species=Dianthus_caryophyllus]:1, beet[species=Beta_vulgaris]:1):1):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( cleome[species=Tarenaya_hassleriana]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1):1, ( ( ( lupine[species=Lupinus_angustifolius]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, (( poplar[species=Populus_trichocarpa]:1, willow[species=Salix_purpurea]:1):1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1)';


    #      '( ((ostreococcus_t[species=Ostreococcus_tauri]:1, ostreococcus_l[species=Ostreococcus_lucimarinus]:1):1, (chlamydomonas[species=Chlamydomonas_reinhardtii]:1, volvox[species=Volvox_carteri]:1):1):1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( (norway_spruce[species=Picea_abies]:1, loblolly_pine[species=Pinus_taeda]:1):1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, (banana[species=Musa_acuminata]:1,( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( nelumbo[species=Nelumbo_nucifera]:1, ( ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, (monkey_flower[species=Mimulus_guttatus]:1, bladderwort[species=Utricularia_gibba]:1):1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, beet[species=Beta_vulgaris]:1):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( cleome[species=Tarenaya_hassleriana]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1):1, ( ( ( lupine[species=Lupinus_angustifolius]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1)';

    # older, lacks cleome, utricularia, nelumbo, etc.
    # '( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, beet[species=Beta_vulgaris]:1):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';
  }
  my $species_tree_parser =
    CXGN::Phylo::Parse_newick->new( $species_tree_newick, 0 );
  my $species_tree =
    $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

  $species_tree->set_missing_species_from_names()
    ;			      # get species from name if species undef
  $species_tree->impose_branch_length_minimum();
  $species_tree->collapse_tree();
  $species_tree->get_root()->recursive_implicit_names();
  $species_tree->get_root()->recursive_implicit_species();

  return $species_tree;
}

sub put_min_branchlengths_into_newick{
  my $newick = shift;
  my $min_bl = shift || 0.00001;
  $newick =~ s/:0[.]0+([,)])/:$min_bl$1/g;
  return $newick;
}


# read in next alignment from STDIN (alfastas file format);
# expect Id ... line, then fasta lines (w gaps),
# then blank line (or EOF?)
# if there is no next alignment (i.e. hit EOF before Id ... line), then returns ('', '', undef, 0)
sub next_align{
 #  my $fh = shift;
   my $id_line = '';
   my $idline_famsize;
   my $count_famsize= 0;
   my $fasta = '';
   while (<>) {
      next if(/^\s*#/);
      if (/^Id.*fam_size:\s*(\d+)/) {
         $id_line = $_;
         $idline_famsize = $1;
         last;
      } else {
         warn "Expected next line to start with 'Id ', but didn't.\n";
      }
   }
   while (<>) {
      last if(/^\s*$/);         # blank line - done
      $count_famsize++ if(/^>/);
      $fasta .= $_;
   }
   return ($id_line, $fasta, $idline_famsize, $count_famsize);
}
