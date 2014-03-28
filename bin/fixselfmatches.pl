#!/usr/bin/perl -w
uses strict;

my $prev_id1 = undef;
while(<>){
my ($id1, $id2, $eval) = split(" ", $_);
if(!defined $prev_id1 or ($id1 ne $prev_id1)){
# new id1. first line should have $id2 eq $id1
  if($id1 ne $id2){
    print "$id1  $id1  0\n";
  }
}
    print $_;
}
