#!/usr/bin/perl -w
use strict;

# find the max e-value in an abc file

my $eval_column = shift || 2; # default is 2 (0-based) for abc


my $maxeval = -1;
while(<>){
my @cols = split(" ", $_);
my $eval = $cols[$eval_column];
	if($eval > $maxeval){
	$maxeval = $eval;
print "$maxeval \n";
}
}
