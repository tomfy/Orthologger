#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $verbose = 0;
GetOptions(
	   'verbose!'      => \$verbose,
	  );

my $fam_count = 0;
my $fam_size = -100;
my ($id, $idline_famsize) = ('no_id', 0);

while (<>) {
   if (/^Id\s+(\S+).*fam_size:\s+(\d+)/) {
      if($id eq 'no_id'){
         print "id               fam_count   fam_size  id_famsize  \n";
      }else{
         my $OK = ($fam_size == $idline_famsize);
         print "$id   $fam_count   $fam_size   $idline_famsize   ",
           $OK? 'OK' : 'XX', "\n"  if($verbose or ($fam_size != $idline_famsize));
      }      
      $id = $1;
      $idline_famsize = $2;
      $fam_size = 0;
      $fam_count++;
   } elsif (/^>/) {
      $fam_size++;
   }
}
print "$id   $fam_count   $fam_size   $idline_famsize   ", 
  ($fam_size == $idline_famsize)? 'OK' : 'XX', " \n" if($verbose or ($fam_size != $idline_famsize));
print "# $fam_count families. Done.\n";
