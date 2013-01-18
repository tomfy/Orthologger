#!/usr/bin/perl -w
use strict;

while(<>){
	my @cols = split(" ", $_);
print $cols[0], "  ", $cols[1], "  ", $cols[10] , "\n";
}
