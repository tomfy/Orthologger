#!/usr/bin/perl -w
use strict;

my $maxeval = shift || 1e-6;

while(<>){
my ($id1, $id2, $eval) = split(" ", $_);
print $_ if($eval <= $maxeval);
}
