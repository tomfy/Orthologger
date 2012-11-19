#!/usr/bin/perl -w
use strict;

no lib '/home/tomfy/bin';
no lib '/home/tomfy/Orthologger/bin';
no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';
use FindBin qw($Bin);
use lib "$Bin/../lib"; 

#use File::Slurp qw ( slurp );

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
my $verbose = shift || 0;

if(!(defined $fasta_filename and defined $cluster_filename and defined $idfile)){
	die "clusters_ids_2_fasta.pl called with no arguments.\n Usage: cluster_ids_2_fasta.pl x.fasta clusterfile idfile. \n x.fasta includes fasta for all sequences in clusterfile, clusterfile has 1 cluster (cluster id & sequence ids) per line. idfile has sequence ids of interest in leftmost col.\n";
}
my %clusterid_seqids = (); # keys are cluster names, values are string of space separated ids
my %seqid_clusterids = (); # keys are ids, values are strings of space separated cluster names (should be just one cluster)

#my @ids = split("\n", slurp($idfile));

open my $fh_id_in, "<$idfile";
my @ids = <$fh_id_in>;
close $fh_id_in;

my %selected_clusters = (); # keys: cluster ids; values: strings with the ids of interest which are in the cluster 
my $no_cluster_count = 0;
for my $id (@ids){ # loop of the ids
$id = ($id =~ /^\s*(\S+)/)? $1: undef;
next if(!defined $id); # handle blank lines in id file
#$id = $1; # just use the first col. 
my $grep_id = "'$id'"; # enclose it in single quotes, or | character will cause a problem in grep command.
#print STDERR "[", $grep_id, "]\n";
#print STDERR "before grep.  id: $id, grep_id:[", $grep_id, "]\n";
# IMGA| -> IMGA_
# $id =~ s/^IMGA[|]/IMGA_/;

	my $cluster = `grep $grep_id $cluster_filename`; # get line (i.e. the cluster) with $id
#print "after grep. \n\n"; #cluster $cluster.\n\n";
	my $cluster_id;
	if($cluster =~ /\S/){ # if non-empty
		$cluster =~ s/([^:]+:)\s*//; # remove initial part with cluster id, 
			$cluster_id = $1; # and put cluster id in scalar $cluster_id
		$cluster_id =~ s/:\s*$//; # remove final colon
	$cluster_id =~ s/[(,]/_/g; $cluster_id =~ s/[)]//;
	$cluster_id =~ s/[\s:]//g; # remove whitespace
	}else{
		$no_cluster_count++;
		$cluster_id = "no_cluster_$no_cluster_count";
		$cluster = $id;
	}
	if(exists $selected_clusters{$cluster_id}){
		my $old_cluster = $selected_clusters{$cluster_id};
		die "cluster inconsistency: old:\n[$old_cluster]\nnew:[$cluster]" if($old_cluster ne $cluster);
	}else{
		$selected_clusters{$cluster_id} = $cluster; # store this cluster id with it string of seq ids.
	}

$clusterid_seqids{$cluster_id} .= "$id "; # concatenate ids (space separated)
$seqid_clusterids{$id} .= "$cluster_id "; 
}

if ($verbose) {
  foreach (keys %clusterid_seqids) {
    print "[$_] [", $clusterid_seqids{$_}, "]\n"
  }
  foreach (keys %seqid_clusterids) {
    print "[$_] [", $seqid_clusterids{$_}, "]\n";
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

        print "$output_fasta_file_name ", 
	  $clusterid_seqids{$cluster_id}, "\n";  #  "  ", `fasta_stats.pl < $output_fasta_file_name`;
      }
# end



