#!/usr/bin/perl -w
use strict;

my %id_eval = ();
while(<>){
my ($id1, $id2, $eval) = split(" ", $_);
$id_eval{$id1} = $eval;
}

my $log10 = log(10);
for(keys %id_eval){
my $ev = $id_eval{$_};
print "$_   $ev   ", ($ev > 0)? log($ev)/$log10 : -10000, "\n";
}

