#!/usr/bin/perl -w
use strict;

my $ggfilename = shift;

my %id_species = ();
if (defined $ggfilename  and  -f $ggfilename) { # read in id-species info
  open my $fhgg, "<", "$ggfilename";
  while (<$fhgg>) {
    my @ids = split(" ", $_);
    my $species = shift @ids;
    $species =~ s/:$//;		# remove final colon
    #	print "species: $species  \n";
    for (@ids) {
      $id_species{$_} = $species;
    }
  }
} else {
  die "No gene-genome association file name supplied. \n";
}

my %desired_species = ('Medicago_truncatula' => 1,
		       'Arabidopsis_thaliana' => 1,
		       'Selaginella_moellendorffii' => 1,
		       'Oryza_sativa' => 1,
		       'Glycine_max' => 1,
		      );


for (<>) {
  my ($id1, $id2, $e_value) = split(" ", $_);
  if (exists $id_species{$id2}) {
    my $id2_species = $id_species{$id2};
    if (exists $desired_species{$id2_species}) {
      print $_;
    }
  }
}

# the end
