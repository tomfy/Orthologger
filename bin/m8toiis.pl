#!/usr/bin/perl -w
use strict;

my $first_line = <>;
my @cols = split(" ", $first_line);
my ($qid, $id2, $ev) = @cols[0, 1, 10];
print "$qid\n";
print "  $id2 $ev\n";
my ($old_qid, $old_id2) = ($qid, $id2);

while ( my $line = <>) {
   @cols = split(" ", $line);
   ($qid, $id2, $ev) = @cols[0,1,10];
   if ($qid ne $old_qid) {
      print "$qid\n";
      print "  $id2 $ev\n";
      $old_qid = $qid;
      $old_id2 = $id2;
   }elsif($id2 ne $old_id2){
      print "  $id2 $ev\n";
      $old_id2 = $id2;
   }
}

