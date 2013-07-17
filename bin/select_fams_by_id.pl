#!/usr/bin/perl 
use strict;

my $ids_filename = shift; # ids are in 1st columns, other cols ignored
my $fastas_filename = shift;

open my $fh_ids, "<", "$ids_filename";

open my $fh_fastas, "<", $fastas_filename;

my %ids = ();
while(<$fh_ids>){
  next if(/^\s*#/);
  if(/^\s*(\S+)/){
    $ids{$1}++;
  }
}

my $print_this_line = 0;
while(<$fh_fastas>){
  if(/^Id (Medtr\S+)/){
    $print_this_line = (exists $ids{$1})? 1 : 0;
  }
  print if($print_this_line);
}
close $fh_ids;
close $fh_fastas;
