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

my $monocots6 = {
		 'Zea_mays'                => 1, # maize
		 'Brachypodium_distachyon' => 1,
		 'Sorghum_bicolor'         => 1,
		 'Oryza_sativa'            => 1, # rice
		 'Musa_acuminata' => 1,		 # banana
		 'Phoenix_dactylifera' => 1,	 # date palm
		};
my $b4 = {
	  'Amborella_trichopoda' => 1,
	  'Picea_abies' => 1,
	  'Selaginella_moellendorffii' => 1,
	  'Physcometrella_patens' => 1,
#	  'Aquilegia_coerulea' => 1,
	 };
my $bother4 = {
	       'Volvox_carteri' => 1,
	       'Chlamydomonas_reinhardtii' => 1,
	       'Ostreococcus_tauri' => 1,
	       'Ostreococcus_lucimarinus' => 1,
	      };

my $pdicots13 = {
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
my %pd13species_present = ();
my %monocot6species_present = ();
my %b4species_present = ();
my %bother4species_present = ();

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
	    scalar keys %monocotspecies_present, "   ",
	      scalar keys %pd13species_present, "   ",
		scalar keys %monocot6species_present, "   ",
		  scalar keys %b4species_present, "   ",
		    scalar keys %bother4species_present, "   ",
		      "\n";
      }
    }
    %pd7species_present = ();
    %monocotspecies_present = ();
    %pd13species_present = ();
    %monocot6species_present = ();
    %b4species_present = ();
    %bother4species_present = ();
    $idmatches = 0;
    $old_id = $id1;
  }
  $idmatches++;
  #	print "id2: $id2 \n";
  my $species2 = $id_species{$id2};
  if (!defined $species2) {
    print "$id1 $id2 \n";
  }
  $pd7species_present{$species2}++ if(exists $pdicots7->{$species2});
  $monocotspecies_present{$species2}++ if(exists $monocots4->{$species2});
  $pd13species_present{$species2}++ if(exists $pdicots13->{$species2});
  $monocot6species_present{$species2}++ if(exists $monocots6->{$species2});
  $b4species_present{$species2}++ if(exists $b4->{$species2});
  $bother4species_present{$species2}++ if(exists $bother4->{$species2});
  
}
print "$old_id  $idmatches   ", 
  scalar keys %pd7species_present, "   ", 
  scalar keys %monocotspecies_present, "   ",
  scalar keys %pd13species_present, "   ",
  scalar keys %monocot6species_present, "   ",
  scalar keys %b4species_present, "   ",
  scalar keys %bother4species_present, "   ",
  "\n";
