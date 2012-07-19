#!/usr/bin/perl -w
use strict;
use lib '/home/tomfy/bin';

#use FastaSelect;
use File::Slurp qw ( slurp );

# given a file with each cluster represented as a line with cluster id
# followed by tab separated sequence ids,
# and a fasta file with ids and sequences for every sequence in clusters
# generate a fasta file for each cluster

# usage:
# clusters_2_fasta.pl fastafile clusterfile
# find all the clusters that the ids in idfile belong to.
# for each such cluster
# print in fasta format the ids and sequences in fastafile
# which are in the cluster.

my $fasta_filename   = shift;
my $cluster_filename = shift;
my $cluster_max_size = 20;
my $cluster_min_size = 5;
my $min_n_taxa = 5;
my $max_n_taxa = 5;

#my $idfile = shift;
my @clusters = split( "\n", slurp($cluster_filename) );    # each line is a cluster

#my %selected_clusters = (); # keys: cluster ids; values: strings with the ids of interest which are in the cluster
my $singleton_count = 0;
my %seqid_clusterid = ();                                  # keys: seq ids; values: ids of corresponding cluster.
for my $cluster (@clusters) {

    my $cluster_id;
    if ( $cluster =~ /\S/ ) {
        $cluster =~ s/([^:]+:)\s*//;                       # remove initial part with cluster id
        $cluster_id = $1;                                  # and put cluster id in scalar $cluster_id

        # now $cluster is tab separated seq ids in cluster; $cluster_id is the id of the cluster.
    } else {
        $singleton_count++;
        $cluster_id = "singleton_$singleton_count";
        #	$cluster = $id;
    }
    $cluster_id =~ s/[(,]/_/g;
    $cluster_id =~ s/[)]//;
    $cluster_id =~ s/[\s:]//g;                             #remove whitespace
    $cluster_id =~ /_(\d+)taxa/;
    my $n_taxa = $1;
    my %cluster_seq_ids = ();                              # keys: ids; values: count;

    my @seq_ids = split( " ", $cluster );                  # list of ids in the cluster.
    next if ( scalar @seq_ids > $cluster_max_size or scalar @seq_ids < $cluster_min_size );
    next if ($n_taxa < $min_n_taxa or $n_taxa > $max_n_taxa);
    foreach (@seq_ids) {
        s/[(].*[)]//g;                                     # remove the stuff in () (i.e. species) for each id
        $cluster_seq_ids{$_}++;
#	print "seq,cluster ids: [$_] [$cluster_id]\n";
        $seqid_clusterid{$_} = $cluster_id;
    }

    # print join(" ", @seq_ids), "\n";
}    # end storing clusters in %seqid_clusterid

#foreach (keys %seqid_clusterid){ print $_, "  ", $seqid_clusterid{$_}, "\n"; }

my %cluster_fasta = ();
die "couldnt open $fasta_filename \n" unless open my $fh, "<$fasta_filename";
my $cluster_id = 'nocluster';
while ( my $line = <$fh> ) {

    if ( $line =~ /^>(\S+)/ ) {    # id line, get id and decide whether to print this one.
        my $seq_id = $1;
#	print "seq id: $seq_id \n";
        $cluster_id = ($seqid_clusterid{$seq_id} || 'nocluster');
#	print "cluster id: $cluster_id\n";
#	exit if (!defined $cluster_id);
    }
    $cluster_fasta{$cluster_id} .= $line unless $cluster_id eq 'nocluster';    # append the line to the appropriate cluster's value.
}
close $fh;
foreach my $cluster_id ( keys %cluster_fasta ) {

    #my $output_fasta_file_name
    #	print $output_fasta_file_name, "  ", `fasta_stats.pl < $output_fasta_file_name`;
    my $output_fasta_file_name = "$cluster_id.fasta";
    open my $fhout, ">$output_fasta_file_name";

    print $fhout $cluster_fasta{$cluster_id};

}

# end

