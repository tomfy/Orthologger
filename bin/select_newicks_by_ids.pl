#!/usr/bin/perl -w
use strict;

# usage select_lines_by_ids.pl idfile < unselected.newicks > selectedfile
# given a file with ids in first col, and
# a newicks file (with multiple newicks in it), format :
# Id Medtr...
# ( -- newick expresion --  )
# 
# i.e. Id line is followed by newick expression on next line
# write out the newicks of the ids which are found in idfile

my $id_file = shift;

my %ids = ();
my @ids_array = ();
my %id_newick = ();
open my $fh, "<", $id_file;
while(<$fh>){
	my @cols = split(" ", $_);
	my $bs_support = (scalar @cols ge 5)? $cols[4] : -1;
	/^\s*(\S+)/;
	$ids{$1} = $bs_support;
	push @ids_array, $1;
#	print "$1 $bs_support \n";
}
#exit;
while(<>){
	my $id = (/^Id\s*(\S+)\s/)? $1 : undef;
	if(defined $id  and  exists $ids{$id}){
		my $string = '';
		$string .= "bs: " . $ids{$id} . "  $_";
		my $newick_line = <>;
		$string .= "$newick_line \n";
		$id_newick{$id} = $string;
	}
}
#my @sids = sort { $ids{$a} <=> $ids{b} } keys %ids;

for(@ids_array){
	print $id_newick{$_};
}