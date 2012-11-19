#!/usr/bin/perl -w
use strict;

# Read in a set of fasta files (as specified by the argument), 
# each of which represents 1 family. 
# Count how many sequences from each of 13 taxa are found
# in each family.
# then summarize; for each taxon: 
# in how many families is it absent
# what is total number of sequences of that taxon in all families.

my $pattern = shift || '*.fasta';
my @ORfiles = split(" ", `ls $pattern`);

my %taxon_taxid = (	'arabidopsis' => '>AT', 
			'brachypodium' => 'Bradi', 
			'castorbean' => '\d{5}[.]m\d{6}',
			'grape' => 'GSVIV',
			'maize' => 'GRMZM|AC\d{6}',
			'medicago' => '>IMGA[|](Medtr|AC|CU)',
			'papaya' => '>evm', # Carica papaya
			'poplar' => 'POPTR',
			'rice' => 'LOC_Os',
			'selaginella' => 'Selmo',
			'sorghum' => '>Sb',
			'soybean' => 'Glyma',
			'tomato' => 'Solyc'
		  );
my %taxon_abrev = (
		   'arabidopsis' => 'ARAth',
		   'brachypodium' => 'BRAdi',
		   'castorbean' => 'RICco',
		   'grape' => 'VITvi',
		   'maize' => 'ZEAma',
		   'medicago' => 'MEDtr',
		   'papaya' => 'CARpa',
		   'poplar' => 'POPtr',
		   'rice' => 'ORYsa',
		   'selaginella' => 'SELmo',
		   'sorghum' => 'SORbi',
		   'soybean' => 'GLYma',
		   'tomato' => 'SOLly'
		  );

my %taxon_nabsfam = (); #key: taxon; value: number of families with this taxon absent.
my @sorted_taxa = sort { $a cmp $b } keys %taxon_taxid;
my $taxa_string = sprintf("       file        genes  taxa ");
foreach (@sorted_taxa) {
  $taxa_string .= sprintf("%5s ", $taxon_abrev{$_});
}
$taxa_string .= "\n";
print $taxa_string;
my %taxon_allclustercount = ();
my %file_resultstring = ();

# my $OR_files = `ls OR*`;
#my @ORfiles = split(" ", `ls $pattern`);

foreach my $orfile (@ORfiles) {
  open my $fhin, "<$orfile";
  my %taxon_count = ();
  while (<$fhin>) {
    if (/^>/) {
      chomp;
      my $id = $_;
      my $id_matches_taxon_in_set = 0;
      foreach my $taxon (keys %taxon_taxid) {
	my $taxid = $taxon_taxid{$taxon};
	if ($id =~ /$taxid/) {
	  $taxon_count{$taxon}++;
	  $taxon_allclustercount{$taxon}++;
	  $id_matches_taxon_in_set++;
	  last;
	}
      }
      die "in $orfile, couldn\'t associate id $id with taxon. \n" if($id_matches_taxon_in_set == 0); 
    }
  }

  my $ngenes_check = 0;

  my ($orf, $ngenes, $ntaxa) = split("_", $orfile);
  $ngenes =~ s/(\d+).*/$1/;
  $ntaxa =~ s/(\d+).*/$1/;
  my $result_string = sprintf("%14s  %5i %5i ", $orf, $ngenes, $ntaxa );
  foreach (@sorted_taxa) {
    if (! exists $taxon_count{$_}) {
      $taxon_nabsfam{$_}++;
    }
    my $count = (exists $taxon_count{$_})? $taxon_count{$_}: 0;
    $result_string .= sprintf("%5i ", $count);
    $ngenes_check += $count;
  }
  $result_string .= "\n";

  $file_resultstring{$orf} = $result_string;
  warn "$orf. ngenes check NG: $ngenes $ngenes_check \n" if($ngenes ne $ngenes_check);
}

my @sorted_files = sort { orthomcl_file_number($a) <=> orthomcl_file_number($b) } keys %file_resultstring;

foreach (@sorted_files) {
  print "  ", $file_resultstring{$_};
}
print $taxa_string;
print "families with taxon absent:   ";
foreach (@sorted_taxa) {
  my $nabsfam = (exists $taxon_nabsfam{$_})? $taxon_nabsfam{$_}: 0;
  printf("%5i ", $nabsfam);
}
print "\n";

my $total_sequences = 0;
print "taxon total, all clusters:    ";
foreach (@sorted_taxa) {
  my $allclustercount = (exists $taxon_allclustercount{$_})? $taxon_allclustercount{$_}: 0;
  printf("%5i ", $allclustercount);
  $total_sequences += $allclustercount;
}
print "\n";
print "total sequences, all taxa, all clusters: $total_sequences \n";
print "nseqs(taxon)/nseqs(all taxa):  ";
foreach (@sorted_taxa) {
  my $allclustercount = (exists $taxon_allclustercount{$_})? $taxon_allclustercount{$_}: 0;
  printf("%1.3f ", $allclustercount/$total_sequences);
}
print "\n";

# subroutines:

sub orthomcl_file_number{
  my $filename = shift;
  if ($filename =~ /(\d+)/) {
    return $1;
  }
  return -1;
}
