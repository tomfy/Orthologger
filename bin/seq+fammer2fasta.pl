#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CXGN::Phylo::CladeSpecifier;

# read in fammer output
# for each id pair (i.e. id1 id2 )
# get the corresponding sequences (fasta)
# usage example:
#  seq+fammer2fasta.pl -gg ../21species.gg  -fammer_file Mv21_10000.fammerout -fasta_in ../new21species-pep.fasta
# output (fasta sequences for each family) would be file:  Mv21_fam_1000.fastas

my $predefined_taxon_groups =
  { # hashref. keys are names of predef taxon groups; values are hashrefs (keys taxa, values 1)
   '4nonangiosperms' => {
			 'Chlamydomonas_reinhardtii'  => 1,
			 'Physcomitrella_patens'      => 1,
			 'Selaginella_moellendorffii' => 1,
			 'Pinus_taeda'                => 1,
			},
   '8monocots' => {
		   'Phoenix_dactylifera'     => 1, # date palm
		   'Setaria_italica'         => 1, # foxtail millet
		   'Triticum_aestivum'       => 1, # wheat
		   'Hordeam_vulgare'         => 1, # barley
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1 # rice
		  },
   '4monocots' => {
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1 # rice
		  },
  '6monocots' => { # These are the monocots in the 36-species analysis, May 2014
		   'Phoenix_dactylifera'     => 1, # date palm
		   'Musa_acuminata' => 1,	   # banana
	      #	   'Setaria_italica'         => 1, # foxtail millet
	#	   'Triticum_aestivum'       => 1, # wheat
      	   #	   'Hordeum_vulgare'         => 1, # barley
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1 # rice
		  },
   '7dicots' => {
		 'Solanum_lycopersicum' => 1, # tomato
		 'Solanum_tuberosum'    => 1, # potato
		 'Vitis_vinifera'       => 1, # grape
		 'Glycine_max'          => 1, # soy
		 'Populus_trichocarpa'  => 1, # poplar
		 'Ricinus_communis'     => 1, # castor
		 'Cucumis_sativus'      => 1  # cucumber
		},
   '8dicots_incl_papaya' => {
			     'Solanum_lycopersicum' => 1, # tomato
			     'Solanum_tuberosum'    => 1, # potato
			     'Vitis_vinifera'       => 1, # grape
			     'Glycine_max'          => 1, # soy
			     'Populus_trichocarpa'  => 1, # poplar
			     'Ricinus_communis'     => 1, # castor
			     'Cucumis_sativus'      => 1, # cucumber
			     'Carica_papaya'        => 1  # papaya
			    },
'13pdicots' => {
		 'Solanum_lycopersicum' => 1, # tomato
		 #	'Solanum_tuberosum'    => 1, # potato
		 'Vitis_vinifera'       => 1,	# grape
		 'Glycine_max'          => 1,	# soy
		 'Populus_trichocarpa'  => 1,	# poplar
		 'Ricinus_communis'     => 1,	# castor
		 'Cucumis_sativus'      => 1,	# cucumber
		 'Aquilegia_coerulea' => 1,	# columbine
		 'Mimulus_guttatus' => 1,	# monkeyflower
		 'Theobroma_cacao' => 1,
		 'Carica_papaya' => 1,
		 'Tarenaya_hassleriana' => 1,
		 'Lupinus_angustifolius' => 1,
		 'Lotus_japonicus' => 1,
		 #	       'Eucalyptus_grandis' => 1,
		 #	       'Manihot_esculenta' => 1,
		},
   '5brassicas' => {
		    Brassica_rapa           => 1, # turnip
		    Arabidopsis_thaliana    => 1,
		    Arabidopsis_lyrata      => 1,
		    Thellungiella_halophila => 1,
		    Capsella_rubella        => 1
		   },
  };
my $taxon_requirements_string = '13pdicots,7 : 6monocots,3'; # '7dicots,6 : 4monocots,3'; # : Selaginella_moellendorffii,1';
my $gg_filename               = undef; # genome-gene association file
my $fammer_file                  = undef; # blast output in abc format
my $input_fasta_filename      = undef;	  # fasta for all sequences
my $max_eval                  = 1e-8;	  # default.
my $max_family_size           = 10000; # default is just a big number, to let families be just whatever is in abc file.

# Process long cl options
GetOptions(
	   'gg_file=s'           => \$gg_filename,
	   'fammer_file=s'          => \$fammer_file,
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
my $min_n_dicots = 0;

########

my $id_sequence_all = store_fasta_sequences($input_fasta_filename);
#print "XXX: ", (exists $id_sequence_all->{'Medtr6g035395.1'})? 'yes' : 'no', "\n";
#exit;
my $gene_genome     = store_gene_genome_association_info($gg_filename);



open my $fh_fammerout, "<", "$fammer_file";
my $previous_id1    = undef;
my $previous_id2    = undef;
my $output_filename = $fammer_file;

$output_filename =~ s/[.](m8|abc)$//;
$output_filename =~ s/_blastout1?$//;
$output_filename .= "_fam.fastas";

open my $fh, ">", "$output_filename";
my ( $fam_size, $fam_string_head, $fam_string_fasta ) = ( 0, '', '' );
my %taxon_count = ();  # keeps track of 
my ($n_core, $fam_size_b, $x, $y, $z, $match_id, $e_value);
my %ids_stored = ();
while ( my $line = <$fh_fammerout> ) {
  next if($line =~ /^\s*#/);  # skip line if first non-whitespace is #
  if ($line =~ /^\s*\[(.*)\]\s*(.*)/) { # first line of family, has core ids in [ ]
    my @core_ids = split(" ", $1);
    my $rest_of_line = $2;
    ($n_core, $fam_size_b, $x, $y, $z) = split(" ", $rest_of_line);
    for (@core_ids) {
      $ids_stored{$_}++;
      $taxon_count{$gene_genome->{$_}}++;
    }
    $fam_size = scalar @core_ids;
    $fam_string_head = "Id [ " . $1 . " ] ";
  } elsif ($line =~ /^\s*(\S+)\s+(\S+)/) {
    ($match_id, $e_value) = ($1, $2);
    $taxon_count{$gene_genome->{$match_id}}++;
    if ( exists $id_sequence_all->{$match_id}) {
      if (!exists $ids_stored{$match_id} ) {
	if ( $fam_size <= $max_family_size ) {
	  $fam_string_fasta .= ">$match_id \n" . $id_sequence_all->{$match_id} . "\n";
	  $fam_size++;
	}
      } else {
	print STDERR "$match_id already stored.\n";
      }
    } else {
      warn "Id [$match_id] not found in $input_fasta_filename.\n";
    }
    
  } elsif ($line =~ /^\s*$/) { # line has only whitespace - end of family; output the family
    my @taxa = sort keys %taxon_count;
    warn "Fam sizes dont agree: [$fam_size] [$fam_size_b].\n" if($fam_size ne $fam_size_b);
    $fam_string_head .= "fam_size: $fam_size  " . join(",", @taxa) . "\n";
    my $taxon_requirement_satisfied = check_taxon_requirements(\@tax_req_objs, \@taxa);
    if ($taxon_requirement_satisfied) {
      print $fh $fam_string_head, $fam_string_fasta, "\n";
 #     print STDERR "Taxa? OK. $fam_string_head \n";
    }else{
#      print STDERR "Taxa? NG. $fam_string_head \n";
    }
    %ids_stored = ();
    %taxon_count = ();
    $fam_size = 0;
    $n_core = 0;
    $fam_string_fasta = '';
    $fam_string_head = '';
for (@tax_req_objs) {
 #  my $the_CS = CXGN::Phylo::CladeSpecifier->new( $_, $predefined_taxon_groups );
# #  print $the_CS->as_string(), "\n";
#   push @tax_req_objs, $the_CS;
  $_->reset();
}
  } else {			# ???
    die "Line has unexpected format: \n", "$line \n";
  }
}


sub store_fasta_sequences {
  my $fasta_filename = shift;
  my %id_sequence    = ();

  open my $fh, "<", "$fasta_filename";
  my ( $id, $sequence ) = ( undef, '' );
  while ( my $line = <$fh> ) {
    if ( $line =~ /^>(\S+)/ ) {
      if ( defined $id ) {
	#  print "[$id] \n";
	$id_sequence{$id} = $sequence;

	#	print "ID: $id \n";

      }
      $id       = $1;
      $sequence = '';
    } else {
      $line =~ s/\s*(\S+)\s*/$1/; # remove whitespace at beginning, end of line.
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
    $genome =~ s/:$//;		# remove the colon.
    for my $the_id (@cols) {
      $gene_genome{$the_id} = $genome;
    }
  }
  return \%gene_genome;
}

sub check_taxon_list { # check list of taxa to see how many monocot species are
  # present, and whether selaginella is present.
  my $taxa = shift;

  my $monocot_count = 0;	# we want to have 3 of these monocots:
  $monocot_count += 1 if ( $taxa =~ /Oryza_sativa/ );
  $monocot_count += 1 if ( $taxa =~ /Sorghum_bicolor/ );
  $monocot_count += 1 if ( $taxa =~ /Brachypodium_distachyon/ );
  $monocot_count += 1 if ( $taxa =~ /Zea_mays/ );

  my $dicot_count = 0;		# we want 6 of these 7 dicots:
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
#      print "$taxon \n";
      
      if ($tax_req->store($taxon)) { # store the taxon and return 1 if requirements have been met.
	$n_satisfied++;
#		print "$n_satisfied \n";
	# print "taxon, n_satisfied: $taxon,  $n_satisfied \n";
	last;
      }
    }
  }
#  print "XXXXX: $n_satisfied, ", scalar @$taxon_requirements, " \n\n\n\n";
  return $n_satisfied == scalar @$taxon_requirements; # true if all the given requirements are satisfied.
}
