#!/usr/bin/perl -w
use strict;

my $min_AM_mono = shift; 
my $min_AM_dicot = shift;
my $min_AM = shift;

while (<>) {
   next if(/^\s*#/);
   my @cols = split(" ", $_);
   my $AM_mono_count = $cols[2];
   my $AM_ambmono_count = $cols[1] + $cols[2];
   my $AM_dicot_count = $cols[3];
   my $AM_ambdicot_count = $cols[1] + $cols[3];
   my $AM_count = $cols[4];

   if 
     # (
     #   ($AM_ambmono_count >= $min_AM_mono) and
     #   ($AM_dicot_count >= $min_AM_dicot) and
     #   ($AM_count >= $min_AM)
     #  ) 
     ( (($AM_ambmono_count >= $min_AM_mono) and
        ($AM_dicot_count >= $min_AM_dicot) ) 
       or
       (($AM_ambdicot_count >= $min_AM_dicot) and
        ($AM_mono_count >= $min_AM_mono) )
     ) {
          print;
       }
}

