#!/usr/bin/perl -w
use strict;

my $max_fam_size = shift || 250;
my %querymatch_count = ();
my %query_matchcount = ();
while(<>){
my @cols = split(" ", $_);
my ($id1, $id2) = @cols[0,1];
$query_matchcount{$id1}++;
next if($query_matchcount{$id1} > $max_fam_size);
$querymatch_count{$id1}++ if($id1 eq $id2);
}
print scalar keys %querymatch_count, "  ", scalar keys %query_matchcount, "\n";
