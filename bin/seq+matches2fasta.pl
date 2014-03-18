#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use CXGN::Phylo::CladeSpecifier;

no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {  # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ;  # the directory containing this script (i.e. Orthologger/bin )
  $libdir = $bindir . '/../lib'; 
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;

use CXGN::Phylo::CladeSpecifier;

# read in blast output
# for each id pair (i.e. id1 id2 )
# get the corresponding sequences (fasta)
# usage example:
#  seq+matches2fasta.pl -gg ../21species.gg  -abc_file Mv21_10000.abc -fasta_in ../new21species-pep.fasta
# output (fasta sequences for each family) would be file:  Mv21_fam_1000.fastas

# my %species_count = (
# 		     # Medicago
# 		     'Medicago_truncatula' => 0,
# 		     # 4 monocots
# 		     'Zea_mays' => 0, 'Brachypodium_distachyon' => 0, 'Oryza_sativa' => 0, 'Sorghum_bicolor' => 0,
# 		     # 8 positive dicots
# 		     'Solanum_lycopersicum' => 0, 'Solanum_tuberosum' => 0, 'Vitis_vinifera' => 0, 'Glycine_max' => 0,
# 		     'Populus_trichocarpa' => 0, 'Ricinus_communis' => 0, 'Cucumis_sativus' => 0, 'Carica_papaya' => 0,
# 		     # 6 negative dicots
# 		     'Arabidopsis_thaliana' => 0, 'Arabidopsis_lyrata' => 0, 'Brassica_rapa' => 0, 'Capsella_rubella' => 0,
# 		     'Thellungiella_halophila' => 0, 'Beta_vulgaris' => 0,
# 		     # Selaginella
# 		     'Selaginella_moellendorffii' => 0, 
# 		     # 3 other
# 		     'Amborella_trichopoda' => 0, 'Chlamydomonas_reinhardtii' => 0, 'Physcomitrella_patens' => 0,
# 		    );
#print "KNOWN SPECIES: ", join(", ", sort keys %species_count), "\n";

my $predefined_taxon_groups =
  {    # hashref. keys are names of predef taxon groups; values are hashrefs (keys taxa, values 1)
    '4nonangiosperms' => {
        'Chlamydomonas_reinhardtii'  => 1,
        'Physcomitrella_patens'      => 1,
        'Selaginella_moellendorffii' => 1,
        'Pinus_taeda'                => 1,
    },
    '8monocots' => {

        'Phoenix_dactylifera'     => 1,    # date palm
        'Setaria_italica'         => 1,    # foxtail millet
        'Triticum_aestivum'       => 1,    # wheat
        'Hordeam_vulgare'         => 1,    # barley
        'Zea_mays'                => 1,    # maize
        'Brachypodium_distachyon' => 1,
        'Sorghum_bicolor'         => 1,
        'Oryza_sativa'            => 1     # rice
    },
    '4monocots' => {
        'Zea_mays'                => 1,    # maize
        'Brachypodium_distachyon' => 1,
        'Sorghum_bicolor'         => 1,
        'Oryza_sativa'            => 1     # rice
    },
    '7dicots' => {
        'Solanum_lycopersicum' => 1,       # tomato
        'Solanum_tuberosum'    => 1,       # potato
        'Vitis_vinifera'       => 1,       # grape
        'Glycine_max'          => 1,       # soy
        'Populus_trichocarpa'  => 1,       # poplar
        'Ricinus_communis'     => 1,       # castor
        'Cucumis_sativus'      => 1        # cucumber
    },
    '8dicots_incl_papaya' => {
        'Solanum_lycopersicum' => 1,       # tomato
        'Solanum_tuberosum'    => 1,       # potato
        'Vitis_vinifera'       => 1,       # grape
        'Glycine_max'          => 1,       # soy
        'Populus_trichocarpa'  => 1,       # poplar
        'Ricinus_communis'     => 1,       # castor
        'Cucumis_sativus'      => 1,       # cucumber
        'Carica_papaya'        => 1        # papaya
    },
    '5brassicas' => {
        Brassica_rapa           => 1,      # turnip
        Arabidopsis_thaliana    => 1,
        Arabidopsis_lyrata      => 1,
        Thellungiella_halophila => 1,
        Capsella_rubella        => 1
    },
  '6negatives' =>  {
        Brassica_rapa           => 1,      # turnip
        Arabidopsis_thaliana    => 1,
        Arabidopsis_lyrata      => 1,
        Thellungiella_halophila => 1,
        Capsella_rubella        => 1,
		  Beta_vulgaris => 1,      # beet
    },
  };
my $default_taxon_requirements_string = '7dicots,6; 4monocots,3'; # ; Selaginella_moellendorffii,1';
my $taxon_requirements_string = $default_taxon_requirements_string;
my $gg_filename               = undef;                                                    # genome-gene association file
my $abc_file                  = undef;                                                    # blast output in abc format
my $input_fasta_filename      = undef;                                                    # fasta for all sequences
my $max_eval                  = 100;                                                     # default is big - 
my $max_family_size           = 10000;   # default is just a big number, to let families be just whatever is in abc file.
my $output_filename = undef;

# Process long cl options
GetOptions(
    'gg_file=s'           => \$gg_filename,
    'abc_file=s'          => \$abc_file,
    'fasta_infile=s'      => \$input_fasta_filename,
    'max_eval=s'          => \$max_eval,
    'max_family_size=i'   => \$max_family_size,
    'taxon_requirement=s' => \$taxon_requirements_string,
    'output_filename=s' => \$output_filename,
);
my $using_default_taxon_requirements = $taxon_requirements_string eq $default_taxon_requirements_string;

print STDERR "OUTPUT FILENAME: $output_filename \n";
my @tax_reqs = split( ":", $taxon_requirements_string );

my @tax_req_objs = ();
for (@tax_reqs) {
#   print "in loop over tax_reqs, $_ \n";
      my $the_CS = CXGN::Phylo::CladeSpecifier->new( $_, $predefined_taxon_groups );
      print $the_CS->as_string(), "\n";
 push @tax_req_objs, $the_CS;
}
my $min_n_monocots = 3;
my $min_n_dicots = 6;

# exit;
########

my $id_sequence_all = store_fasta_sequences($input_fasta_filename);
my $gene_genome     = store_gene_genome_association_info($gg_filename);

open my $fh_blast, "<", "$abc_file";
my $previous_id1    = undef;
my $previous_id2    = undef;
if(!defined $output_filename){ # if no output filename specified on CL, make a file name
$output_filename = $abc_file; # from $abc_file by replacing .abc with .fastas
$output_filename =~ s/[.](m8|abc)$/.fastas/;
}
#print STDERR "output filename:  $output_filename \n";
#exit;
open my $fh, ">", "$output_filename";
my ( $fam_size, $fam_string_head, $fam_string_fasta ) = ( 0, '', '' );
my %taxon_count = ();
while ( my $line = <$fh_blast> ) {
    my @cols = split( " ", $line );
    my ( $id1, $id2, $eval ) = @cols[ 0, 1, 2 ];
    next if ( $eval > $max_eval );
    next if ( defined $previous_id1 and ( $id1 eq $previous_id1 and $id2 eq $previous_id2 ) );

    if ( ( !defined $previous_id1 ) or ( $id1 ne $previous_id1 ) ) {
        my @taxa = sort keys %taxon_count;
        my $cs_taxa = join( ",", @taxa );
#	count_species(\@taxa, \%species_count);
        $fam_string_head .= "fam_size: $fam_size  $cs_taxa\n";
        my ( $n_dicots, $n_monocots, $selaginella_present ) = check_taxon_list($cs_taxa);
	my $old_OK = ($n_dicots >= $min_n_dicots and $n_monocots >= $min_n_monocots); #  and !$selaginella_present);

	my $taxon_requirement_satisfied = check_taxon_requirements(\@tax_req_objs, \@taxa);

	if ($using_default_taxon_requirements and ($old_OK ne $taxon_requirement_satisfied)) {
	  warn "Old, new taxon requirements satisfied:  [$old_OK]  [$taxon_requirement_satisfied] n_dicots: $n_dicots n_monocots: $n_monocots \n";
	}

        # print "XXX: $cs_taxa    [$n_monocots]   [$selaginella_present].\n";
        if ( defined $previous_id1 ) {
	  if ($taxon_requirement_satisfied){ # $n_dicots >= $min_n_dicots and $n_monocots >= 3 and !$selaginella_present ){
            print $fh "$fam_string_head";
            print $fh "$fam_string_fasta";
            print $fh "\n";
	  }
        }

        # $fam_string_head .= "fam_size: $fam_size " . join(",", @taxa) . "\n";
        # #    print "$previous_id1  $fam_size  ", scalar @taxa, "\n";
        # print $fh "$fam_string_head", "$fam_string_fasta", "\n" if(defined $previous_id1);

        #    print $fh "Id $id1 family: \n";
        %taxon_count      = ();
        $fam_string_head  = "Id $id1 family. ";
        $fam_string_fasta = '';
        $previous_id1     = $id1;
        $fam_size         = 0;
    }

    if ( exists $id_sequence_all->{$id2} ) {
        if ( $fam_size <= $max_family_size ) {

            #    print $fh
            $fam_string_fasta .= ">$id2 \n" . $id_sequence_all->{$id2} . "\n";
            $fam_size++;

            #    $previous_id1 = $id1;
            $previous_id2 = $id2;
            my $taxon = $gene_genome->{$id2};

            #  print "$id2  $taxon \n";
            $taxon_count{$taxon}++;
        }
    }
    else {
        warn "Id $id2 not found in $input_fasta_filename.\n";
    }
  } # loop over lines of blast output (abc)
my @taxa = sort keys %taxon_count;
my $cs_taxa = join( ",", @taxa );
$fam_string_head .= "fam_size: $fam_size  $cs_taxa\n";
my ( $n_dicots, $n_monocots, $selaginella_present ) = check_taxon_list($cs_taxa);
my $old_OK = ($n_dicots >= $min_n_dicots and $n_monocots >= $min_n_monocots); #  and !$selaginella_present);

my $taxon_requirement_satisfied = check_taxon_requirements(\@tax_req_objs, \@taxa);
	if($using_default_taxon_requirements and ($old_OK ne $taxon_requirement_satisfied)){
warn "Old, new taxon requirements satisfied:  [$old_OK]  [$taxon_requirement_satisfied] n_dicots: $n_dicots n_monocots: $n_monocots \n";
}
#print " [$old_OK]  [$taxon_requirement_satisfied]  $n_dicots $n_monocots \n";


  # print "XXX: $cs_taxa    [$n_monocots]   [$selaginella_present].\n";
  if ( defined $previous_id1 ) {
    if ($taxon_requirement_satisfied){ #  $n_dicots >= $min_n_dicots and $n_monocots >= 3 and !$selaginella_present ){
    print $fh "$fam_string_head";
    print $fh "$fam_string_fasta";
    print $fh "\n";
  }
  }

# my $output_fasta_filename = "$cluster_id.fasta";
# open my $fhout, ">", "$output_fasta_filename";
# die "couldnt open $input_fasta_filename \n" unless open my $fh_in, "<", "$fasta_filename";
# my $print_on = 0;
# while (my $line = <$fh_in>) {
#   if ($line =~ /^>(\S+)/) { # id line, get id and decide whether to print this one.
#     my $id = $1;
#     $print_on = (exists $cluster_seq_ids{$id})? 1: 0;
#   }
#   print $fhout $line if($print_on);
# }
# close $fhout;
# close $fh;
# print "$output_fasta_file_name ",
#   $clusterid_seqids{$cluster_id}, "\n"; #  "  ", `fasta_stats.pl < $output_fasta_file_name`;


sub store_fasta_sequences {
    my $fasta_filename = shift;
    my %id_sequence    = ();

    open my $fh, "<", "$fasta_filename";
    my ( $id, $sequence ) = ( undef, '' );
    while ( my $line = <$fh> ) {
        if ( $line =~ /^>(\S+)/ ) {
            if ( defined $id ) {
                $id_sequence{$id} = $sequence;

                #	print "ID: $id \n";

            }
            $id       = $1;
            $sequence = '';
        }
        else {
            $line =~ s/\s*(\S+)\s*/$1/;    # remove whitespace at beginning, end of line.
            $sequence .= $line;
        }
    }
    if ( defined $id ) {
        $id_sequence{$id} = $sequence;

        #    print "ID: $id \n";
    }
    return \%id_sequence;
}

sub store_gene_genome_association_info {
    my $gg_filename = shift;
    my %gene_genome = ();

    open my $fh, "<", "$gg_filename";
    while (<$fh>) {
        my @cols = split( " ", $_ );
        my $genome = shift @cols;
        $genome =~ s/:$//;    # remove the colon.
        for my $the_id (@cols) {
            $gene_genome{$the_id} = $genome;
        }
    }
    return \%gene_genome;
}

sub check_taxon_list {        # check list of taxa to see
  # how many of the specified monocot species are present,
  # how many of the specified dicot species are present,
  # and whether selaginella is present.
    my $taxa = shift;

    my $monocot_count = 0;    # counts how many of these 4 monocot species are present:
    $monocot_count += 1 if ( $taxa =~ /Oryza_sativa/ );
    $monocot_count += 1 if ( $taxa =~ /Sorghum_bicolor/ );
    $monocot_count += 1 if ( $taxa =~ /Brachypodium_distachyon/ );
    $monocot_count += 1 if ( $taxa =~ /Zea_mays/ );

    my $dicot_count = 0;      # counts how many of these 7 dicot species are present:
    $dicot_count += 1 if ( $taxa =~ /Populus_trichocarpa/ );
    $dicot_count += 1 if ( $taxa =~ /Ricinus_communis/ );
    $dicot_count += 1 if ( $taxa =~ /Solanum_lycopersicum/ );
    $dicot_count += 1 if ( $taxa =~ /Solanum_tuberosum/ );
    $dicot_count += 1 if ( $taxa =~ /Vitis_vinifera/ );
    $dicot_count += 1 if ( $taxa =~ /Glycine_max/ );
    $dicot_count += 1 if ( $taxa =~ /Cucumis_sativus/ );
    
# print "taxon list string: $taxa; m,d counts: $monocot_count, $dicot_count \n";
    return ( $dicot_count, $monocot_count, $taxa =~ (/Selaginella/) ? 1 : 0 );
}

sub check_taxon_requirements {
  my $taxon_requirements = shift; # this is an arrayref of CladeSpecifier objects
  my $taxa               = shift; # array ref
  my $n_satisfied        = 0;
  for my $tax_req (@$taxon_requirements) { # $tax_req is a CladeSpecifier obj.
    $tax_req->reset();
    for my $taxon (@$taxa) {
      if ($tax_req->store($taxon)){ # store the taxon and return 1 if requirements have been met.
	$n_satisfied++;
	last;
      }
    }
    $tax_req->reset();
  }
  return $n_satisfied == scalar @$taxon_requirements; # true if all the given requirements are satisfied.
}

# sub count_species{
#   my $taxa = shift; # array ref
#   my $species_count = shift; # hashref
#   for(@$taxa){
#     if(exists $species_count->{$_}){
#       $species_count->{$_}++;
#     }else{
#       print STDERR "Species  [$_]  unknown.\n";
#     }
#   }
#   my @specs = sort keys %$species_count;
#   for(@specs){
# printf("%4i ", $species_count->{$_});
# }printf("\n");
#   for(@specs){ $species_count->{$_} = 0; }
# return 1;
# }
