#!/usr/bin/perl -w
use strict;

my $gg_filename = shift;
my $max_of_one_species = shift || 10;

# hashref. keys: seq ids, values: taxa
my $gene_genome = store_gene_genome_association_info($gg_filename);
my %taxon_count = ();

my $prev_id1 = 'X_X_X_X_X';

while (<>) {
  my ($id1, $id2, $eval) = split(" ", $_);
  my $taxon;
  if ($id1 ne $prev_id1) {
    my @staxa = sort keys %taxon_count;
    printf("%18s ", $prev_id1);
    for (@staxa) {
      printf("%4i ", $taxon_count{$_});
    }
    print "\n";
    for (keys %taxon_count) {
      $taxon_count{$_} = 0;
    }
    # %taxon_count = ();
    $prev_id1 = $id1;
    #  print "$taxon  $id1 \n";
  }
  $taxon = $gene_genome->{$id2};
  $taxon_count{$taxon}++;
  if ($taxon_count{$taxon} <= $max_of_one_species) {
    #   print $_;
  }
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
