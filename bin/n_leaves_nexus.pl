#!/usr/bin/perl -w
use strict;

my @nexi = `ls *.nexus`;
#print join(", ", @nexi), "\n";


for my $nexus_filename (@nexi){
	my $n_leaves = 0;
	my %sp_count = ();
	chomp $nexus_filename;
#	print "[$nexus_filename]\n";
	open my $fh, "<", "$nexus_filename" or die "couldnt open $nexus_filename for reading. \n";
	while(<$fh>){
		$n_leaves++ if(/color=/);
		$sp_count{$1}++ if(/__(\S+)'\[/);
	}
	close $fh;
	my $qid = $nexus_filename;
	$qid =~ s/[.]nexus//;
	print "$qid  $n_leaves ", scalar keys %sp_count, "\n";
}
