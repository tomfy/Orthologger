#!/usr/bin/perl -w
use strict;

my %id1_count = ();

while(<>){
my @cols = split(" ", $_);
my $id1 = $cols[0];
$id1_count{$id1}++
}

print "number of ids: ", scalar keys %id1_count, "\n";
