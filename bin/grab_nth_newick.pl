#!/usr/bin/perl -w
use strict;

my $n = shift || 1;
my $qid;
my $fam_size;
my $count = 0;
while(<>){
if(/^Id\s+(\S+)\s+family[.]\s+fam_size:\s+(\S+)/){
$qid = $1;
$fam_size = $2;
}
if(/^ML\s*(.*)/){
	$count++;
	if($count == $n){
	print STDERR "$qid  $fam_size\n";
	print $1;
}
}
}
