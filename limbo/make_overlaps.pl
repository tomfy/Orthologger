#!/usr/bin/perl 
use strict;
use Getopt::Std;
use List::Util qw ( min max sum );

use lib '/home/tomfy/Orthologger/lib';
use Overlap;
# use Orthologger;
# use TlyUtil qw ( order_newick );

use lib '/home/tomfy/cxgn/cxgn-corelibs/lib';
use CXGN::Phylo::File;
use CXGN::Phylo::Parser;

my $pattern = shift || '*.fasta';
my @files = `ls $pattern`;
my $nongap_fraction = shift || 0.8; 

foreach my $input_file (@files){

my $output_file = $input_file;
$output_file =~ s/[.]{2}\///g;
print "input file: $input_file, output_file: $output_file \n";
#### Get the alignment:
my $align_string =  `cat $input_file`;
# fixes to $align_string:
$align_string =~ s/IMGA[|]/IMGA_/g; #pipe in id is a problem; '|' -> '_'.
$align_string =~ s/(>[^|]+)[|][^\n]*\n/$1\n/g; # delete from first pipe to end of line.
my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this.
$align_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg; 
# print "align string: $align_string\n";

# construct an overlap object.
my $overlap_obj = Overlap->new($align_string, $nongap_fraction);
my $overlap_fasta_string = $overlap_obj->overlap_fasta_string('');

open my $fh, ">$output_file";
print $fh $overlap_fasta_string;
close $fh;

}
