#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CXGN::Phylo::CladeSpecifier;

# read in blast output
# for each id pair (i.e. id1 id2 )
# get the corresponding sequences (fasta)
# usage example:
#  seq+matches2fasta.pl -gg ../21species.gg  -abc_file Mv21_10000.abc -fasta_in ../new21species-pep.fasta
# output (fasta sequences for each family) would be file:  Mv21_fam_1000.fastas

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
my $taxon_requirements_string = '7dicots,6 : 4monocots,3'; # : Selaginella_moellendorffii,1';
my $gg_filename               = undef;                                                    # genome-gene association file
my $abc_file                  = undef;                                                    # blast output in abc format
my $input_fasta_filename      = undef;                                                    # fasta for all sequences
my $max_eval                  = 1e-8;                                                     # default.
my $max_family_size           = 10000;   # default is just a big number, to let families be just whatever is in abc file.

# Process long cl options
GetOptions(
    'gg_file=s'           => \$gg_filename,
    'abc_file=s'          => \$abc_file,
    'fasta_infile=s'      => \$input_fasta_filename,
    'max_eval=s'          => \$max_eval,
    'max_family_size=i'   => \$max_family_size,
    'taxon_requirement=s' => \$taxon_requirements_string
);

my @tax_reqs = split( ":", $taxon_requirements_string );

my @tax_req_objs = ();
for (@tax_reqs) {
   
      my $the_CS = CXGN::Phylo::CladeSpecifier->new( $_, $predefined_taxon_groups );
      print $the_CS->as_string(), "\n";
 push @tax_req_objs, $the_CS;
}
my $min_n_dicots = 6;

########

my $id_sequence_all = store_fasta_sequences($input_fasta_filename);
my $gene_genome     = store_gene_genome_association_info($gg_filename);

open my $fh_blast, "<", "$abc_file";
my $previous_id1    = undef;
my $previous_id2    = undef;
my $output_filename = $abc_file;

# print STDERR $abc_file, "\n";
$output_filename =~ s/[.](m8|abc)$//;
$output_filename =~ s/_blastout1?$//;
$output_filename .= "_fam.fastas";

# print STDERR "$output_filename \n";
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
        $fam_string_head .= "fam_size: $fam_size  $cs_taxa\n";
        my ( $n_dicots, $n_monocots, $selaginella_present ) = check_taxon_list($cs_taxa);
my $taxon_requirement_satisfied = check_taxon_requirements(\@tax_req_objs, \@taxa);
my $old_OK = ($n_dicots >= $min_n_dicots and $n_monocots >= 3 and !$selaginella_present);
#print " [$old_OK]  [$taxon_requirement_satisfied] \n";

        # print "XXX: $cs_taxa    [$n_monocots]   [$selaginella_present].\n";
        if ( defined $previous_id1 ) {
            print $fh "$fam_string_head";
            print $fh "$fam_string_fasta"
              if ( $n_dicots >= $min_n_dicots and $n_monocots >= 3 and !$selaginella_present );
            print $fh "\n";
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
my $taxon_requirement_satisfied = check_taxon_requirements(\@tax_req_objs, \@taxa);
my $old_OK = ($n_dicots >= $min_n_dicots and $n_monocots >= 3 and !$selaginella_present);
#print " [$old_OK]  [$taxon_requirement_satisfied] \n";


  # print "XXX: $cs_taxa    [$n_monocots]   [$selaginella_present].\n";
  if ( defined $previous_id1 ) {
    print $fh "$fam_string_head";
    print $fh "$fam_string_fasta" if ( $n_dicots >= $min_n_dicots and $n_monocots >= 3 and !$selaginella_present );
    print $fh "\n";
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

sub check_taxon_list {        # check list of taxa to see how many monocot species are
                              # present, and whether selaginella is present.
    my $taxa = shift;

    my $monocot_count = 0;    # we want to have 3 of these monocots:
    $monocot_count += 1 if ( $taxa =~ /Oryza_sativa/ );
    $monocot_count += 1 if ( $taxa =~ /Sorghum_bicolor/ );
    $monocot_count += 1 if ( $taxa =~ /Brachypodium_distachyon/ );
    $monocot_count += 1 if ( $taxa =~ /Zea_mays/ );

    my $dicot_count = 0;      # we want 6 of these 7 dicots:
    $dicot_count += 1 if ( $taxa =~ /Populus_trichocarpa/ );
    $dicot_count += 1 if ( $taxa =~ /Ricinus_communis/ );
    $dicot_count += 1 if ( $taxa =~ /Solanum_lycopersicum/ );
    $dicot_count += 1 if ( $taxa =~ /Solanum_tuberosum/ );
    $dicot_count += 1 if ( $taxa =~ /Vitis_vinifera/ );
    $dicot_count += 1 if ( $taxa =~ /Glycine_max/ );
    $dicot_count += 1 if ( $taxa =~ /Cucumis_sativus/ );

    return ( $dicot_count, $monocot_count, $taxa =~ (/Selaginella/) ? 1 : 0 );
}

sub check_taxon_requirements {
  my $taxon_requirements = shift; # this is an arrayref of CladeSpecifier objects
  my $taxa               = shift; # array ref
  my $n_satisfied        = 0;
  for my $tax_req (@$taxon_requirements) {
    for my $taxon (@$taxa) {
      if ($tax_req->store($taxon)){ # store the taxon and return 1 if requirements have been met.

	$n_satisfied++;
	# print "taxon, n_satisfied: $taxon,  $n_satisfied \n";
	last;
      }
    }
  }
  return $n_satisfied == scalar @$taxon_requirements; # true if all the given requirements are satisfied.
}
