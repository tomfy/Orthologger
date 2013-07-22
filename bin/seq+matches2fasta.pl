#!/usr/bin/perl -w
use strict;

# read in blast output
# for each id pair (i.e. id1 id2 )
# get the corresponding sequences (fasta)

my $blast_out_file = shift;	# id1 id2 .... eval (abc format)
my $input_fasta_filename = shift;
my $gg_filename = shift;
my $max_eval = shift || 0.000001;
my $max_family_size = shift || 200;

my $min_n_dicots = 0;

my $id_sequence_all = store_fasta_sequences($input_fasta_filename);
my $gene_genome = store_gene_genome_association_info($gg_filename);

open my $fh_blast, "<", "$blast_out_file";
my $previous_id1 = undef;
my $previous_id2 = undef;
my $output_filename = $blast_out_file;
# print STDERR $blast_out_file, "\n";
$output_filename =~ s/[.](m8|abc)$//;
$output_filename =~ s/_blastout1?$//;
$output_filename .=  "_fam.fastas";
# print STDERR "$output_filename \n";
#exit;
open my $fh, ">", "$output_filename";
my ($fam_size, $fam_string_head, $fam_string_fasta) = (0, '', '');
my %taxon_count = ();
while (my $line = <$fh_blast>) {
  my @cols = split(" ", $line);
  my ($id1, $id2, $eval) = @cols[0,1,2];
  next if($eval > $max_eval);
  next if(defined $previous_id1  and  ($id1 eq $previous_id1  and   $id2 eq $previous_id2));

  if ((! defined $previous_id1)  or ($id1 ne $previous_id1) ) {
    my @taxa = sort keys %taxon_count;
my $cs_taxa = join(",", @taxa);
$fam_string_head .= "fam_size: $fam_size  $cs_taxa\n";
my ($n_dicots, $n_monocots, $selaginella_present) = check_taxon_list($cs_taxa);
# print "XXX: $cs_taxa    [$n_monocots]   [$selaginella_present].\n";
if (defined $previous_id1) {
  print $fh "$fam_string_head";
  print $fh "$fam_string_fasta" if($n_dicots >= $min_n_dicots  and  $n_monocots >= 3  and  $selaginella_present);
  print $fh "\n";
}


    # $fam_string_head .= "fam_size: $fam_size " . join(",", @taxa) . "\n";
    # #    print "$previous_id1  $fam_size  ", scalar @taxa, "\n";
    # print $fh "$fam_string_head", "$fam_string_fasta", "\n" if(defined $previous_id1);

    #    print $fh "Id $id1 family: \n";
    %taxon_count = ();
    $fam_string_head = "Id $id1 family. ";
    $fam_string_fasta = '';
    $previous_id1 = $id1;
    $fam_size = 0;
  }

  if (exists $id_sequence_all->{$id2}) {
    if ($fam_size <= $max_family_size ) {

      #    print $fh 
      $fam_string_fasta .= ">$id2 \n" . $id_sequence_all->{$id2} . "\n";
      $fam_size++;
      #    $previous_id1 = $id1;
      $previous_id2 = $id2;
      my $taxon = $gene_genome->{$id2};
      #  print "$id2  $taxon \n";
      $taxon_count{$taxon}++;
    }
  } else {
    warn "Id $id2 not found in $input_fasta_filename.\n";
  }
}
my @taxa = sort keys %taxon_count;
my $cs_taxa = join(",", @taxa);
$fam_string_head .= "fam_size: $fam_size  $cs_taxa\n";
my ($n_dicots, $n_monocots, $selaginella_present) = check_taxon_list($cs_taxa);
# print "XXX: $cs_taxa    [$n_monocots]   [$selaginella_present].\n";
if (defined $previous_id1) {
  print $fh "$fam_string_head";
  print $fh "$fam_string_fasta" if($n_dicots >= $min_n_dicots  and  $n_monocots >= 3 and $selaginella_present);
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

sub store_fasta_sequences{
  my $fasta_filename = shift;
  my %id_sequence = ();

  open my $fh, "<", "$fasta_filename";
  my ($id, $sequence) = (undef, '');
  while (my $line = <$fh>) {
    if ($line =~ /^>(\S+)/) {
      if (defined $id) {
	$id_sequence{$id} = $sequence;
	#	print "ID: $id \n";

      }
      $id = $1;
      $sequence = '';
    } else {
      $line =~ s/\s*(\S+)\s*/$1/; # remove whitespace at beginning, end of line.
      $sequence .= $line;
    }
  }
  if (defined $id) {
    $id_sequence{$id} = $sequence;
    #    print "ID: $id \n";
  }
  return \%id_sequence;
}

sub store_gene_genome_association_info{
  my $gg_filename = shift;
  my %gene_genome = ();

  open my $fh, "<", "$gg_filename";
  while (<$fh>) {
    my @cols = split(" ", $_);
    my $genome = shift @cols;
    $genome =~ s/:$//;		# remove the colon.
    for my $the_id (@cols) {
      $gene_genome{$the_id} = $genome;
    }
  }
  return \%gene_genome;
}

sub check_taxon_list{ # check list of taxa to see how many monocot species are
  # present, and whether selaginella is present.
  my $taxa = shift;

  my $monocot_count = 0;	# we want to have 3 of these monocots:
  $monocot_count += 1 if($taxa =~ /Oryza_sativa/);
  $monocot_count += 1 if($taxa =~ /Sorghum_bicolor/);
  $monocot_count += 1 if($taxa =~ /Brachypodium_distachyon/);
  $monocot_count += 1 if($taxa =~ /Zea_mays/);

  my $dicot_count = 0;		# we want 6 of these 7 dicots:
  $dicot_count += 1 if($taxa =~ /Populus_trichocarpa/);
  $dicot_count += 1 if($taxa =~ /Ricinus_communis/);
  $dicot_count += 1 if($taxa =~ /Solanum_lycopersicum/);
  $dicot_count += 1 if($taxa =~ /Solanum_tuberosum/);
  $dicot_count += 1 if($taxa =~ /Vitis_vinifera/);
  $dicot_count += 1 if($taxa =~ /Glycine_max/);
  $dicot_count += 1 if($taxa =~ /Cucumis_sativus/);

  return ($dicot_count, $monocot_count, $taxa =~ (/Selaginella/)? 1 : 0);
}
