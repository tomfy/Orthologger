#!/usr/bin/perl 
use strict;

# 1st argument: an id list file; ids are in leftmost column.
# 2nd argument: a m8 format blastout file (or abc)

my $ids_filename = shift; # ids are in 1st columns, other cols ignored
my $m8_filename = shift;

open my $fh_ids, "<", "$ids_filename" or die "couldn't open $ids_filename for reading.\n";

open my $fh_m8, "<", "$m8_filename" or die "couldn't open $m8_filename for reading.\n";

my %ids = ();
while(<$fh_ids>){
  next if(/^\s*#/);
  if(/^\s*(\S+)/){
    $ids{$1}++;
  }
}
close $fh_ids;

my $print_this_line = 0;
while(<$fh_m8>){
  if(/^\s*(\S+)/){
    $print_this_line = (exists $ids{$1})? 1 : 0;
  }
  print if($print_this_line);
}
close $fh_m8;
