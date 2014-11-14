#!/usr/bin/perl -w
use strict;

my %id_count = ();
my $gg_filename = shift;

my $gene_genome     = store_gene_genome_association_info($gg_filename);

my %species_idcount = ();     # keys: species, values id:count hashref

my $id;
while (<>) {
  if (/^>(\S+)/) {
    #  my ($id1, $id2, $eval) = split(" ", $_);
    $id = $1;
  } elsif (/^(\S+)\s+(\S+)\s+\S+\s*$/) { # abc line
    $id = $2;
  } else {
    next;
  }
  my $species = $gene_genome->{$id};
  $id_count{$id}++;
  if (exists $species_idcount{$species}) {
    $species_idcount{$species}->{$id}++;
  } else {
    $species_idcount{$species} = {$id => 1};
  } 
  #  $id_count{$id2}++;
}

  print "number of ids in families: ", scalar keys %id_count, "\n";
while (my ($sp, $idc) = each %species_idcount) {
  print "$sp  ", scalar keys %$idc, "\n";
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
