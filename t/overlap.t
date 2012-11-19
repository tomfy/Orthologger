#!/usr/bin/perl -w
use strict;
use warnings FATAL => 'all';

# tests for Overlap Module
use Test::More tests=> 2;


use File::Basename 'dirname';
use Cwd 'abs_path';
my ($bindir, $libdir);
BEGIN{
$bindir = dirname(abs_path(__FILE__)); # this has to go in Begin block so happens at compile time
$libdir = $bindir . '/../lib';
$libdir = abs_path($libdir); # collapses the bin/../lib to just lib
}
use lib $libdir;

use CXGN::Phylo::Parser;
use CXGN::Phylo::Overlap;

my $file = `ls OR*`;
$file =~ s/\s*$//;
print "file: [$file]\n";
my $overlap_obj = CXGN::Phylo::Overlap->new($file);

ok( defined $overlap_obj, 'new() returned something.' );
isa_ok( $overlap_obj, 'CXGN::Phylo::Overlap' );

