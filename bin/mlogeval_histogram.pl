#!/usr/bin/perl -w
use strict;

my $big_positive_value = 300;
my $log10 = log(10.0);
while(<>){
my ($id1, $id2, $e_val) = split(" ", $_);
my $mlogeval = ($e_val >= 1e-300 )? -1*log($e_val)/$log10 : $big_positive_value;

print "$id1 $id2 $mlogeval \n";
} 
