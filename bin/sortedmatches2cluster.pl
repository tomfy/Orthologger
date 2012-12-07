#!/usr/bin/perl -w
use strict;

# gets up to $size_limit best matches read in from stdin (assumed to
# already be sorted with earlier lines being better matches (smaller e-value) )
# and write out (to stdout) in cluster file format, i.e. family_name: id1 id2 (etc.)
#i.e. all on one line, family name is followed by colon, and then by whitespace-separated
# sequence ids of sequences in family.

my $size_limit = shift || 100;

my %id_count = (); 
my $line_count = 0;
my ($id1, $id2, $eval);
my $family_name = "fam_";
my $ids_string = "";
while(<>){
($id1, $id2, $eval) = split(" ", $_);
$id_count{$id1}++; $line_count++;
if($id_count{$id1} != $line_count){
	warn "Ids in lefthand col not all the same? At line $line_count.\n";
}
if($id_count{$id1} == 1){
	$family_name .= "$id1";
#	$ids_string .= "$id1 ";
}
$ids_string .= "$id2 ";
last if($line_count >= $size_limit);
}
print "$family_name:  $ids_string \n";


