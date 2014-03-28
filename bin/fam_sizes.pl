#!/usr/bin/perl -w
use strict;

my $fam_count = 0;
my $fam_size = -100;
while(<>){
	if(/^Id M/){
	print "$fam_count   $fam_size \n";
	$fam_size = 0;
	$fam_count++;
}elsif(/^>/){
	$fam_size++;
}
}
print "$fam_count   $fam_size \n";
