#!/usr/bin/perl -w
use strict;

my %id1_matchcount = ();
my %qid_selfmatch = ();
my $qmatch_count = 0;
while(<>){
my @cols = split(" ", $_);
my $id1 = shift @cols;
my $id2 = shift @cols;
if($id1 eq $id2){
  $qid_selfmatch{$id1} = 1;
}
$id1_matchcount{$id1}++;
}

print "Families: ", scalar keys %id1_matchcount, "; ",
  "queries match selves: ", 
  scalar keys %qid_selfmatch, " \n";


for(keys %id1_matchcount){
  if(!exists $qid_selfmatch{$_}){
    print "query $_ doesn't match self. ??\n";
  }
}
