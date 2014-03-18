#!/usr/bin/perl -w
use strict;

# get 1 newicks expression from a newicks file
# read from stdin, write to stdout
# usage example: get1newick.pl 'Medtr1g007170.1'  < p1.newicks

my $id = shift;
my $v = shift || 0;
my $njml = shift || 'ML';

my $id_found = 0;

while (<>) {
  last if(/^Id\s+\Q$id\E/);
}
while (<>) {
  if (/^\Q$njml\E\s+(.*)/) {
    print "# id: $id, tree algorithm: $njml \n";
    print "$1 \n";
    exit;
  }
}
die("Found no  $njml  tree with id $id in it.\n");

