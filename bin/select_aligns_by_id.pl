#!/usr/bin/perl -w
use strict;

my $idfile = shift;
my $alf_filename_pattern = shift;
open my $fhid, "<", $idfile || die "couldnt open $idfile for reading.\n";

my %id_count = ();
for (<$fhid>) {
   next if(/^\s*#/);
   if (/^\s*(\S+)/) {
      my $id = $1;
      $id_count{$id}++;
   }
}
print "$alf_filename_pattern \n";

my @files = `ls $alf_filename_pattern`;
print STDERR join("; ", @files), "\n";
exit;
for my $alf_file (@files) {
   $alf_file =~ s/\s+$//;
   open my $fh_alf, "<", $alf_file or die "couldnt open $alf_file for reading. \n";
   #my $idpattern = shift;
   my $print = 0;
   while (<$fh_alf>) {
      if (/^Id\s+(\S+)/) {
         if (exists $id_count{$1}) {
            $print = 1;
            #			print "$_ turning on print \n";
         } else {          # not one of desired ids, turn off printing
            $print = 0;
         }
         #		$print = 0;
      }
      #else{
      print $_ if($print);
      #	}
   }
}

