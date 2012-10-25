#!/usr/bin/perl -w
use strict;

my $cluster_file = shift;
my $fasta_file = shift;

# write out all the clusters which have an id found in the $fasta_file

my %clusterid_idstring = (); # keys cluster ids; values string of sequence ids.
my %seqid_clusterid = ();
my %selected_clusters = ();
open my $fh, "<$cluster_file";
while(<$fh>){
	s/^(OR[^:]+:)\s*//;
	my $cluster_id = $1;
	s/\s+$//; # delete final whitespace 
		$clusterid_idstring{$cluster_id} = $_;
	my @ids = split(" ", $_);

	foreach my $id (@ids){
	$id =~ s/[(].*[)]//;
		$seqid_clusterid{$id} = $cluster_id;
	}
}
close $fh;

my ($n_found, $n_not_found) = (0, 0);

open $fh, "<$fasta_file";
while(<$fh>){
	my $seq_id;
	if(/^>(\S+)/){
		$seq_id = $1;

		if(exists $seqid_clusterid{$seq_id}){
			my $cluster_id = $seqid_clusterid{$seq_id};
			$selected_clusters{$cluster_id}++;
	$n_found++;
		}else{
	$n_not_found++;
		#	warn "sequence id: $seq_id not found in cluster file $cluster_file\n";
		}
	}
}
foreach my $cluster (keys %selected_clusters){
	print "$cluster  ", $clusterid_idstring{$cluster}, "\n";
}
exit;

foreach (keys %seqid_clusterid){
if(/m00/){
print $_, "  ", $seqid_clusterid{$_}, "\n";
}
}
print "found/not found: $n_found, $n_not_found \n";
