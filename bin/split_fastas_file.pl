#!/usr/bin/perl -w
use strict;

my $fasta_file = shift;

my $n_parts = shift || 2; # number of pieces to split into
my $wcout = `grep '^>' $fasta_file | wc`;
#print $wcout;
#my $wcout = `wc $grepout`;
#print $wcout, "\n";
#exit;
my @cols = split(" ", $wcout); # `grep '^>' $fasta_file  |  wc`)
	
my $n_seqs = $cols[0];
open my $fh, "<", "$fasta_file";

my $n_sequences_in_each_part = int($n_seqs/$n_parts) + 1;
my $target_sequence_count = $n_sequences_in_each_part;
print "number of sequences: $n_seqs ; number of parts: $n_parts number of sequences per part: $target_sequence_count \n";

my $i_part = 1;
my $output_filename = $fasta_file . ".part$i_part";
open my $fh_out, ">", $output_filename;
my $sequence_count = 0; 
while(<$fh>){
	if(/^>/){
	$sequence_count++;
#	print "$sequence_count  $_\n";
	}elsif(/^Id /){
	if($sequence_count > $target_sequence_count){
	close $fh_out;
	$i_part++;
	$output_filename = $fasta_file . ".part$i_part";
	$target_sequence_count += $n_sequences_in_each_part;
	open $fh_out, ">", $output_filename;
 }
}
	print $fh_out  $_;
}
close $fh_out;
close $fh;
