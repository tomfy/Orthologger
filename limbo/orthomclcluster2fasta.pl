#!/usr/bin/perl -w
use strict;
use lib '/home/tomfy/bin';
use FastaSelect;

# use to get the cluster containing a particular sequence id,
# and output fasta for that cluster.
# Orthomcl gives a file with the ids in each cluster.
# Here we give it a fasta file which also has the sequences.
# This finds the cluster with argument 'id' in it, gets all the ids
# in that cluster, and outputs the fasta for these ids (assuming the ids
# are present in the fastafile.)

# usage:
# orthomclcluster2fasta.pl fastafile clusterfile id
# print in fasta format the ids and sequences in fastafile 
# which are in the cluster in clusterfile which contains id

my $fasta_filename = shift;
my $cluster_filename = shift;
my $id = shift;
my %seqids = ();

my $cluster = `grep $id $cluster_filename`;
# print "$cluster\n";

$cluster =~ s/([^:]+:)\s*//;
my $cluster_id = $1;
my @seq_ids = split(" ", $cluster);
foreach (@seq_ids){
	s/[(].*[)]//;

	$seqids{$_}++; 
}
# print join(" ", @seq_ids), "\n";

open my $fh, "<$fasta_filename";
my $print_on = 0;
while(<$fh>){
	if(/^>(\S+)/){ # id line, get id and decide whether to print this one.
	my $id = $1;
	$print_on = (exists $seqids{$id})? 1: 0;
}
print if($print_on);
}
	

	
