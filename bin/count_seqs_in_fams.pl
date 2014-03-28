#!/usr/bin/perl -w
use strict;

#

my %id_count = ();
while(<>){
my ($id1, $id2, $eval) = split(" ", $_);
$id_count{$id1}++;
$id_count{$id2}++;
}

print "number of ids in families: ", scalar keys %id_count, "\n";
