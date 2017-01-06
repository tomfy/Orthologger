#!/usr/bin/perl -w
use strict;
use SpeciesCriteria;

my $infile = shift;

my $species_str = 'Medicago_truncatula,';
$species_str .= 'Phaseolus_vulgaris, Solanum_lycopersicum,';
$species_str .= 'Oryza_sativa, ';
#$species_str .= 'Brachypodium_distachyon,';
#$species_str .= 'Amborella_trichopoda,';
$species_str .= 'Capsella_rubella,';
#$species_str .= 'Beta_vulgaris,';
$species_str =~ s/\s+//g;
my %species_present = ();
for my $sp (split(",", $species_str)) {
   $species_present{$sp}++;
}

my $SpCritObj = SpeciesCriteria->new($infile);
my $OK = $SpCritObj->check_criteria(keys %species_present);
print "OK: ", $OK? '1' : '0', "\n";
