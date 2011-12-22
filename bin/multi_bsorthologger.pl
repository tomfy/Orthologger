#!/usr/bin/perl -w
use strict;

my $ls_pattern = shift;
my $n_bootstrap = shift || 400;
my $q_species = shift || 'castorbean';
my $min_support = shift || 0.15;
my $seed = shift || 1123456;
my $files_string = `ls $ls_pattern`;
my @filenames = split(" ", $files_string);

print "filenames: ", join("\n", @filenames), "\n";

my $count = 0;
foreach my $align_filename (@filenames){
  print "$count  $align_filename \n";
  my $output_filename = $align_filename;
$output_filename =~ s/^[.]{2}\///; # remove ../ if present
  $output_filename =~ s/fasta$//;
  $output_filename .= "bs$n_bootstrap" . 'ortholog';
  my $cl = "perl ~/MHarrisonProject/bin/bootstrap_ortholog.pl -i $align_filename";
$cl .= " -T ML -k -n 1 -N $n_bootstrap -S $seed -r mindl -m $min_support -q $q_species > $output_filename";
 system "$cl";
  $seed += 100;
  $count++;
}
