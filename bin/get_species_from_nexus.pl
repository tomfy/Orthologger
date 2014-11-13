#!/usr/bin/perl -w
use strict;

my %species = ();
my %sequence_species = ();

while(<>){
last if(/^\s*taxlabels/);
}

while(<>){
  if(/^\s*['](\S+)__(\w+_\w+)[']/){
my ($seq_id, $sp) = ($1, $2);
    $species{$sp}++;
$sequence_species{$seq_id} = $sp;
  }
}

# my @ssp = sort {$species{$b} <=> $species{$a}}  keys %species;

# for(@ssp){
#   print "$_  ", $species{$_}, "\n";
# }

while(my ($id, $sp) = each %sequence_species){
  printf("%40s  %30s  %8i\n", $id, $sp, $species{$sp});
}
