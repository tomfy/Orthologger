#!/usr/bin/perl -w
use strict;

my $nmin_pd7 = shift;
my $nmin_mono4 = shift || 0;
my $nmin_pd13 = shift || 0;
my $nmin_mono6 = shift || 0;
my $nmin_b4 = shift || 0;
my $nmin_bother4 = shift || 0;
my $nmin_b8 = shift || 0;

while(<>){
my ($id, $famsize, $npd7, $nmono4, $npd13, $nmono6, $nb4, $nbother4) = split(" ", $_);
if($npd7 >= $nmin_pd7 and $nmono4 >= $nmin_mono4 and $npd13 >= $nmin_pd13 and $nmono6 >= $nmin_mono6
	and $nb4 >= $nmin_b4 and $nbother4 >= $nmin_bother4
and ($nb4 + $nbother4)>= $nmin_b8){

print $_;
}
}
