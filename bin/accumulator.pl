#!/usr/bin/perl -w
use strict;

# take a file where each line has several numbers representing
# e.g. numbers of events of certain kinds in some time interval,
# and generate a corresponding file with numbers which are
# sums up to this point.

# usage:  accumulator.pl < infile

my $decay_length = shift || undef;

my $accumulate_col_0 = shift || 0;
#print "decay length: $decay_length \n";

my @accum = split(" ", <>);

while (<>) {
  my @cols = split(" ", $_);
  foreach my $i (0..scalar @cols-1) {

    if ($i == 0 and !$accumulate_col_0) {
      $accum[$i] = $cols[$i]; # col 0 is likely to be the number of steps or some such - don't want to accumulate
    } else {
      if (defined $decay_length) {
	$accum[$i] = (1.0 - 1.0/$decay_length) * $accum[$i] + $cols[$i]/$decay_length; 
      } else {
	$accum[$i] += $cols[$i];
      }
    }
  }
  print join(" ", @accum), "\n";
}

