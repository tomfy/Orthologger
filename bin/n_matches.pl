#!/usr/bin/perl -w
use strict;

my $min_monocots = 3;
my $min_pdicots = 6;

my $monocots4 = {
		 'Zea_mays'                => 1, # maize
		 'Brachypodium_distachyon' => 1,
		 'Sorghum_bicolor'         => 1,
		 'Oryza_sativa'            => 1	# rice
		};
my $pdicots7 = {
		'Solanum_lycopersicum' => 1, # tomato
		'Solanum_tuberosum'    => 1, # potato
		'Vitis_vinifera'       => 1, # grape
		'Glycine_max'          => 1, # soy
		'Populus_trichocarpa'  => 1, # poplar
		'Ricinus_communis'     => 1, # castor
		'Cucumis_sativus'      => 1  # cucumber
	       };

my $ggfilename = shift;
open my $fhgg, "<", "$ggfilename";
my %id_species = ();
while (<$fhgg>) {
  my @ids = split(" ", $_);
  my $species = shift @ids;
  $species =~ s/:$//;		# remove final colon
  #	print "species: $species  \n";
  for (@ids) {
    $id_species{$_} = $species;
  }
}
# exit;

my $idmatches = 0;
my %pd7species_present = ();
my %monocotspecies_present = ();
my $dummy_id = 'Not a reasonable id';
my $old_id = $dummy_id;
while (<>) {			# read in blast results, abc format.
  my ($id1, $id2, $e_value) = split(" ", $_);
  if ($id1 ne $old_id) {
    if ($old_id ne $dummy_id) {
      my $npdicots = scalar keys %pd7species_present;
      my $nmonocots = scalar keys %monocotspecies_present;
      if (1  or  ($npdicots >= $min_pdicots  and  $nmonocots >= $min_monocots)) {
	print "$old_id  $idmatches   ", 
	  scalar keys %pd7species_present, "   ", 
	    scalar keys %monocotspecies_present, "\n";
      }
    }
    %pd7species_present = ();
    %monocotspecies_present = ();
    $idmatches = 0;
    $old_id = $id1;
  }
  $idmatches++;
  #	print "id2: $id2 \n";
  my $species2 = $id_species{$id2};
  if(!defined $species2){
    print "$id1 $id2 \n";
  }
  $pd7species_present{$species2}++ if(exists $pdicots7->{$species2});
  $monocotspecies_present{$species2}++ if(exists $monocots4->{$species2});
}
	print "$old_id  $idmatches   ", 
	  scalar keys %pd7species_present, "   ", 
	    scalar keys %monocotspecies_present, "\n";
