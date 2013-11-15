#!/usr/bin/perl -w
use strict;

# take out from select_from_cladesout.pl
# keep only IDs for which ML is good (i.e. '1'),
# output the ID and number of accepted bootstraps.

my %id_nbsg = ();
while(<>){
	my @cols = split(" ", $_);
	next if($cols[2] ne '1');
	my ($id, $n_bs_good) = ($cols[0], $cols[4]);
	$id_nbsg{$id} = $n_bs_good;
}

my @sorted_ids = sort {$id_nbsg{$b} <=> $id_nbsg{$a}} keys %id_nbsg;

for(@sorted_ids){
	print "$_ ", $id_nbsg{$_}, "\n";
}
