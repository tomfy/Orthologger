#!/usr/bin/perl -w
use strict;

my %id_count = ();
my %qid_count = ();
while (<>) {
  my ($id1, $id2, $eval) = split(" ", $_);
  $id_count{$id1}++;
  $id_count{$id2}++;
  $qid_count{$id1}++;
}

print "There are ", scalar keys %qid_count, " families (queries) and ", scalar keys %id_count, " total matches.\n";
