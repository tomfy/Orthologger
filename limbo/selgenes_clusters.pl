#!/usr/bin/perl -w
use strict;

#  usage: selgenes_clusters.pl < selids
# selids is file with 1 id per line
# run it from 

my %cluster_selid = (); # key: cluster file name, value: selected ids in that cluster.
while(my $selid = <>){
	$selid =~ s/\s//g; # delete whitespace 	
		my $cluster = `grep $selid OR*`;
	if(! $cluster =~ /\S/){
		$cluster = `grep $selid single*`;
	}
	if(! $cluster =~ /\S/){
		warn "id: $selid not found in OR* of single* files!?!?!?\n";
	}else{
	$cluster =~ s/\s+$//;
		$cluster =~ /^([^.]+)/;
		my $cluster_file = $1;
		$cluster_selid{$cluster_file} .= "$selid ";

	}
}

my @skis = sort { clustersize($b) <=> clustersize($a) } keys %cluster_selid;

foreach (@skis){
	print $_, "  ", $cluster_selid{$_}, "\n";

}


sub clustersize{
	my $id = shift;
	if($id =~ /_(\d+)genes/){
		return $1;
	}else{
		return 1;

	}
}
