#!/usr/bin/perl -w
use strict;

# simple script to run FastTree
# and write output to files with name based on input file name

my $pattern = shift || "*.fasta";
my @files = `ls $pattern`;

foreach my $input_file (@files){
	chomp $input_file;
# my $input_file = shift;
	my $output_file = $input_file;
	$output_file =~ s/[.]{2}\///g;

	my $FT_cl = "FastTree -wag -gamma -bionj $input_file 2> $output_file.out";
	print $FT_cl, "\n";

	my $fasttree_newick_out = `$FT_cl`;

	open my $fh, ">$output_file.newick";
	print $fh "$fasttree_newick_out";
	close $fh;

}
