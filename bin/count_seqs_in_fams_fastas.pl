#!/usr/bin/perl -w
use strict;

my %id_count = ();
while (<>) {
  if(/^>(\S+)/){
#  my ($id1, $id2, $eval) = split(" ", $_);
    my $id = $1;
 $id_count{$id}++;
#  $id_count{$id2}++;
}
}

print "number of ids in families: ", scalar keys %id_count, "\n";
