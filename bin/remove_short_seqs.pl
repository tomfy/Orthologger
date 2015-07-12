#!/usr/bin/perl -w
use strict;

my $min_seqlength = shift || 32;
my $seqlength = 0;
my $sequence = '';
my $idline = <>;
my $id = '';
if ($idline =~ /^>(\S)/) {
   $id = $1;
}
my ($count, $short_count) = (0, 0);
# print STDERR "XXX: $idline";
print STDERR "min seqlength: $min_seqlength \n";
while (<>) {
   if (/^>(\S+)/) {
      $count++;
      if ($seqlength >= $min_seqlength) {
         print "$idline";
         print "$sequence";
      } else {
         # print STDERR "$id, sequence is too short: $seqlength \n";
         $short_count++;
      }
      if ($count % 1000 == 0) {
         print STDERR "$count  $short_count \r";
      }
      $seqlength = 0;	
      $sequence = '';
      $idline = $_;
      $id = $1;
   } else {                     # sequence line
      my $the_line = $_;
      #	print "A: ", length $the_line, " [$the_line] \n";
      $the_line =~ s/[*\s]+$//; # remove final * and whitespaces
      #	print "B: ", length $the_line, " [$the_line] \n";
      #	print "Line length: ", length $the_line, "\n";
      $seqlength += length $the_line;
      $sequence .= $_;
   }

}
print STDERR "all sequences: $count  short sequences: $short_count \n";

print $idline;
print $sequence;
