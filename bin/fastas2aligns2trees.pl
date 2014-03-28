#!/usr/bin/perl -w
use strict;

my $multifamily_fasta_file = shift;

open my $fh, "<", "$multifamily_fasta_file";

my $fasta_string = '';
while(<$fh>){
  if(/^Id (\S+) family/){
    my $id_line = $_;
    $fasta_string = '';
  }elsif(/^\s*$/){
    
