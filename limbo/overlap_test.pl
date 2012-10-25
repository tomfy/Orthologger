#!/usr/bin/perl -w
use strict;
use lib '/home/tomfy/MHarrisonProject/lib';
use Overlap;

my $infile = shift;
my $fraction = shift || 0.8;
my $overlap_obj = Overlap->new($infile, $fraction);

#print "fasta string for input alignment: \n", $overlap_obj->align_fasta_string(), "\n";

for (1..2){
print $overlap_obj->bootstrap_overlap_fasta_string();
print "\n";
}
