#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $verbose = 0;
GetOptions(
	   'verbose!'      => \$verbose,
	  );

my ($id, $fs, $numb_eq);

while (<>) {
   if (/^FT /) { 
      $numb_eq = $_ =~ tr/=/=/; # each leaf has [species= , so count '='s to count leaves.
   } elsif (/^Id (\S+)\s.*fam_size:\s+(\S+)/) {
      $id = $1; $fs = $2;
   } elsif (/^\s*$/) {
      my $OK = ($fs == $numb_eq);
      print "$id  $fs  $numb_eq ", 
        $OK? "OK" : "XX", "\n" if($verbose or ($fs != $numb_eq));
   }
}
