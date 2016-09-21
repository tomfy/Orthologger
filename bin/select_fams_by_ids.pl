#!/usr/bin/perl 
use strict;

# 1st argument: an id list file; ids are in leftmost column.
# 2nd argument: a fastas format file with (possibly aligned) families each in fasta format. 
# [ e.g.:
# Id Medtr...  # the query sequence id which this family is based on.
# >At...       # the corresponding family fasta
# ABDSFDASDFA...
# >Another_id
# AAGDFDSfDSDF...
#                      # end of this family
# Id Medtr...        # next family query seq. id.
# etc. end of example. ]
# output the fasta for all families with query id in the id list file.

my $ids_filename = shift; # ids are in 1st columns, other cols ignored
my $fastas_filename = shift;

open my $fh_ids, "<", "$ids_filename";

open my $fh_fastas, "<", $fastas_filename  or  die "couldnt open $fastas_filename for reading.\n";

my %ids = ();
while(<$fh_ids>){
  next if(/^\s*#/);
  if(/^\s*(\S+)/){
    $ids{$1}++;
  }
}

my $print_this_line = 0;
while(<$fh_fastas>){
  if(/^Id\s+(\S+)/){
   #  print "$1 \n";
    $print_this_line = (exists $ids{$1})? 1 : 0;
  }
  print if($print_this_line);
}
close $fh_ids;
close $fh_fastas;
