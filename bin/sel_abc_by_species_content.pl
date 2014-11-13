#!/usr/bin/perl -w
use strict;

my $gg_filename;
my $abc_filename;

GetOptions(
	   'gg_file=s'           => \$gg_filename,
	   'abc_file=s'          => \$abc_filename,
	  );

my $gg = store_gene_genome_association_info($gg_filename);

my @required_species = ('Setaria_italica', 'Sorghum_bicolor', 'Zea_mays', 'Lupinus_angustifolius');

my $line = <>;
my @cols = split (" ", $line);
my ($id1, $id2, $ev) = @cols[0..2];
my $species = $gg->{$id2};
my $old_id1 = $id1;
my %species_count = ();
$species_count{$species}++;
my $fam_string = "$id1  $id2  $ev\n";
while ($line = <>) {
  @cols = split (" ", $line);
  ($id1, $id2, $ev) = @cols[0..2];
  if ($id1 ne $old_id1) {	# print the old id
    my $required_species_ok = 1;
    for (@required_species) {
      $required_species_ok = 0 if(!exists $species_count{$_});
    }
    print "$fam_string" if($required_species_ok);
    $old_id1 = $id1;
    $fam_string = '';
  }

  $fam_string = "$id1  $id2  $ev\n";
  

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
