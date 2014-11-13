#!/usr/bin/perl -w
use strict;

# arg is file with desired ids in first col
# sel_abc_by_ids.pl  idfilename < abcfile.abc
# print out all lines of abcfile.file which have id (query id) in $idfile
my $idfile = shift;
open my $fhid, "<$idfile";

# store ids in has
my %id_count = ();
for(<$fhid>){
next if(/^\s*#/);
if(/^\s*(\S+)/){
	$id_count{$1}++;
}
}

while(<>){ # read in line of abc file
next if(/^\s*#/);
my ($id1, $id2, $eval) = split(" ", $_);
if(exists $id_count{$id1}){
	print;
}
}
