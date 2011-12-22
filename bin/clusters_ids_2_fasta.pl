#!/usr/bin/perl -w
use strict;
use lib '/home/tomfy/bin';
use FastaSelect;
use File::Slurp qw ( slurp );

# use to get the clusters containing each of a set of sequence ids,
# and output fasta for each cluster.

# usage:
# clusters_ids_2_fasta.pl fastafile clusterfile idfile
# find all the clusters that the ids in idfile belong to.
# for each such cluster
# print in fasta format the ids and sequences in fastafile 
# which are in the cluster.

my $fasta_filename = shift;
my $cluster_filename = shift;
my $idfile = shift;
my @ids = split(" ", slurp($idfile));
my %selected_clusters = (); # keys: cluster ids; values: strings with the ids of interest which are in the cluster 
my $singleton_count = 0;
for my $id (@ids){ # loop of the ids

	my $cluster = `grep $id $cluster_filename`; # get line (i.e. the cluster) with $id
	my $cluster_id;
	if($cluster =~ /\S/){
		$cluster =~ s/([^:]+:)\s*//; # remove initial part with cluster id
			$cluster_id = $1; # and put cluster id in scalar $cluster_id
	}else{
		$singleton_count++;
		$cluster_id = "singleton_$singleton_count";
		$cluster = $id;
	}
	if(exists $selected_clusters{$cluster_id}){
		my $old_cluster = $selected_clusters{$cluster_id};
		die "cluster inconsistency: old:\n[$old_cluster]\nnew:[$cluster]" if($old_cluster ne $cluster);
	}else{
		$selected_clusters{$cluster_id} = $cluster; # store this cluster id with it string of seq ids.
	}
}

for my $cluster_id (keys %selected_clusters){
	my $cluster = $selected_clusters{$cluster_id};

	my %cluster_seq_ids = ();
	my @seq_ids = split(" ", $cluster); # list of ids in the cluster.
		foreach (@seq_ids){
			s/[(].*[)]//g; # remove the stuff in () (i.e. species) for each id
				$cluster_seq_ids{$_}++;
		}
# print join(" ", @seq_ids), "\n";
	$cluster_id =~ s/[(,]/_/g; $cluster_id =~ s/[)]//;
	$cluster_id =~ s/[\s:]//g; #remove whitespace
		my $output_fasta_file_name = "$cluster_id.fasta";
	open my $fhout, ">$output_fasta_file_name";
	die "couldnt open $fasta_filename \n" unless open my $fh, "<$fasta_filename";
	my $print_on = 0;
	while(my $line = <$fh>){
		if($line =~ /^>(\S+)/){ # id line, get id and decide whether to print this one.
			my $id = $1;
			$print_on = (exists $cluster_seq_ids{$id})? 1: 0;
		}	
		print $fhout $line if($print_on);
	}
	close $fhout;
	close $fh;
	print $output_fasta_file_name, "  ", `fasta_stats.pl < $output_fasta_file_name`;
}
# end



