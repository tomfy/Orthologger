#!/usr/bin/perl -w
use strict;


my %id_count = ();
my $the_col = shift || 1;
my $verbose = shift || undef;
while(<>){
#	if(/^>(\S+)/)
	my @cols = split(" ", $_);
if(scalar @cols >= $the_col){
		$id_count{$cols[$the_col-1]}++;
	}
}
print "#Number of distinct ids:  ", scalar keys %id_count, "\n";
if($verbose){
print join("\n", sort keys %id_count), "\n";
}
