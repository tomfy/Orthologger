#!/usr/bin/perl -w
use strict;
use warnings;

# usage: check_ids_seqs.pl  subset.fasta  superset.fasta
# for each sequence in f1, check that id is also found in  f2.fasta
# and that the sequences agree.
my $subset_file = shift;

my $superset_file = shift; # must be in form in which each sequence is all on one line.
if(!defined $subset_file  or !defined $superset_file or !-f $subset_file or !-f $superset_file){
die "check_ids_seqs.pl requires two fasta files as arguments.\n check_ids_seqs.pl f1.fasta f2.fasta\n";
}

open my $FH_sub, "<$subset_file";
while(<$FH_sub>){

	if(/^>(\S+)/){  # an id line of f1
		my $print_output = '';
		my $idline = $_;
		my $id = $1;
		my $sequence = <$FH_sub>;
		$sequence = uc $sequence;
		$sequence =~ s/[*]*\s*$//; # remove terminal asterisks
#	print "id: $id \n";	
	my $grep_pattern = "'" . $id . "'";
#	print "grep_pattern  $grep_pattern \n";
			my $grepout = ` grep -A 1 $grep_pattern $superset_file `; # check for id in f1
		if($grepout =~ /\S/){ # non-whitespace present in grep output
			my ($line1, $line2) = split("\n", $grepout);
			$line1 =~ /[>]?(\S+)/;
			$line1 = $1;
			$line2 = uc $line2;
			$line2 =~ s/[*]*\s*$//; # remove terminal asterisks
				my $seq_agree = ($sequence eq $line2)? 'YES' : 'NO';
			$print_output = "$id  $line1; sequences agree? $seq_agree\n";
			if($seq_agree eq 'NO'){ print STDERR  "$print_output $sequence\n $line2\n"; }
		print "diff [" . $sequence ^ $line2 . "]\n";		
}else{
			$print_output = "$id ; id not found\n";
		}
		print $print_output;
	}

}
