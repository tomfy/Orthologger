#!/usr/bin/perl -w
use strict;



  my $muscle_input_filename = shift;
  # if ($do_taxonify) {
  #   my $taxonify_cl = $bindir . '/' . 
  #     "taxonify_fasta.pl $fasta_file $gg_filename";
  #   $muscle_input_filename = `$taxonify_cl`;
  # }

  #   $muscle_input_filename =~ s/\s*$//;                # remove terminal whitespace

  #    print "rmuscle input file name: [$muscle_input_filename] \n";
  if ( $muscle_input_filename =~ /^\s*(\S+).*$/) {
    $muscle_input_filename = $1;
  } else {
    die "Filename for input to muscle is all whitespace: [$muscle_input_filename].\n";
  }
  my $align_filename = $muscle_input_filename;
  $align_filename =~ s/[.]fasta/_align.fasta/;
  my $muscle_cl  = "muscle -in $muscle_input_filename -out $align_filename";
  my $muscle_out = `$muscle_cl 2>&1`;
 my $count_warnings = 0;
my @muscle_out_lines = split("\n", $muscle_out);
  for(@muscle_out_lines){
    if(/^[*]{3}\sWARNING\s[*]{3}/){
      $count_warnings++;
      next;
    }
    if(/\S/){
      print;
    }
  }
#print "muscle stdout output: $muscle_out \n";
  print "Done with alignment using muscle. alignment filename: [$align_filename]. $count_warnings warnings. \n";
