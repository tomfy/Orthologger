#!/usr/bin/perl -w
use strict;

# usage select_newicks_by_ids.pl idfile < unselected.newicks > selected.newicks
# given a file with ids in first col, and
# a newicks file (with multiple newicks in it), format :
# Id Medtr...
# ( -- newick expression --  )
# 
# i.e. Id line is followed by newick expression on next line
# write out the newicks of the ids which are found in idfile

my $id_file = shift;

my %ids = ();
my @ids_array = ();
my %id_newick = ();
open my $fh, "<", $id_file;
while(<$fh>){
  next if(/^\s*#/); # ignore comment lines
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
		my $string = $_;
	#	$string .= # "bs: " . $ids{$id} . "$_";
		my $newick_line = <>;
		$string .= "$newick_line \n";
		$id_newick{$id} = $string;
	}
}
#my @sids = sort { $ids{$a} <=> $ids{b} } keys %ids;

for(@ids_array){
  if(exists $id_newick{$_}){
	print $id_newick{$_};
      }else{
	print STDERR "id $_ not a key of id_newick hash.\n";
      }
}
