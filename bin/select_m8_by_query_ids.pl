#!/usr/bin/perl -w
use strict;

#usage: select_m8_by_query_ids.pl xxx_ids < zzz.m8

my $id_filename = shift; # ids in lefmost col, possibly prec3eded by whitepsace.

my %id_count = ();

open my $fh, "<", $id_filename or die "couldn't open $id_filename for reading.\n";
while (<$fh>) {
  if (/^\s*(\S+)/) {
    my $id = $1;
    $id_count{$id}++;
  }
}

while (<>) {
  my @cols = split(" ", $_);
  if (exists $id_count{$cols[0]}) {
    print;
  }
}
