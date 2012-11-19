#!/usr/bin/perl -w
use strict;

# align files with muscle
# give it a pattern to match, by default just *.fasta

my $pattern = shift || '*.fasta';

my @files = `ls $pattern`;


foreach my $infile (@files){
	next if ($infile =~ /align/);
	$infile =~ s/^\s*(\S+)\s*$/$1/;	
	my $outfile = $infile;
	$outfile =~ s/[.]fasta/_align.fasta/;
#	print "infile: $infile, outfile: $outfile \n";
	my $muscle_out = `muscle -in $infile -out $outfile`;

	print $outfile, "\n";
}


