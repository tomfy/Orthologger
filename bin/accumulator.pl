#!/usr/bin/perl -w
use strict;

#my @accum = split(" ", <>);

my $decay_length = shift || undef;
#print "decay length: $decay_length \n";

my @accum = split(" ", <>);

while(<>){
	my @cols = split(" ", $_);
#$accum[0] = $cols[0]; # 0th col not cumulative (ngen)
	foreach my $i (0..scalar @cols-1){
# $accum[$i] += $cols[$i];
		if(defined $decay_length){
			$accum[$i] = (1.0 - 1.0/$decay_length) * $accum[$i] + $cols[$i]/$decay_length; 
		}else{
			$accum[$i] += $cols[$i];
		}
	}
	print join(" ", @accum), "\n";
}
