#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw (min max sum);

no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib'; 
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;

use Phyml;
use CXGN::Phylo::Overlap;
#use CXGN::Phylo::Parser;
#use CXGN::Phylo::BasicTree;
#use CXGN::Phylo::File;
#use CXGN::Phylo::Species_name_map;

my $nongap_fraction = shift || 0.8;

my $alignment_id = 'anonymous';
my $fasta_string = '';
while (<>) {
   if (/^Id\s+(\S+)/) {
      $alignment_id = $1;
   } elsif (/^\s*$/) {
      
      print "Id $alignment_id \n";
      $fasta_string .= "\n";
      my $overlap_obj = CXGN::Phylo::Overlap->new( $fasta_string, $nongap_fraction );
      my $overlap_fasta_string = $overlap_obj->overlap_fasta_string('');
      print "$overlap_fasta_string \n";
      $fasta_string = '';
   } else {
      if (/^>/) {
         $fasta_string .= "\n" . $_;
      } else {
         $_ =~ s/\s+//g; 
         $fasta_string .= $_;
      }
   }
}
# if no whitespace line after fasta, print the overlap string
if ($fasta_string ne '') {
   $fasta_string .= "\n";
   my $overlap_obj = CXGN::Phylo::Overlap->new( $fasta_string, $nongap_fraction );
   my $overlap_fasta_string = $overlap_obj->overlap_fasta_string('');
   print "$overlap_fasta_string \n"; 
}

