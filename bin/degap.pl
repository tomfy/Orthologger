#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
use IPC::Run3;
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
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

my $nongap_fraction = shift || 0.8;
my $rng_seed = 1234567;
my @lines = <>;
my $fasta = join("", @lines);
my $overlap_obj =
  CXGN::Phylo::Overlap->new( $fasta, $nongap_fraction,
			     $rng_seed );
my $overlap_fasta_string =
  $overlap_obj->overlap_fasta_string('');
print $overlap_fasta_string, "\n";
