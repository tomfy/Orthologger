#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $query_id = undef;
my $n = undef;
my $method = 'FT';     # FastTree by default, others are 'NJ', 'PHYML'
my $verbose = 0;
GetOptions(
	   'id|qid=s'           => \$query_id,
	   'n=i' => \$n,
	   'method=s' => \$method,
	   'verbose!' => \$verbose,
	  );
# get 1 newicks expression from a newicks file
# read from stdin, write to stdout
# usage example: get1newick.pl 'Medtr1g007170.1'  < p1.newicks

#my $id = shift;
#my $v = shift || 0;
#my $njml = shift || 'FT';

my $id_found = 0;
my $n_trees_read = 0;

if (defined $query_id) {
  while (<>) {
    last if(/^Id\s+\Q$query_id\E/);
  }
  while (<>) {
    if (/^\Q$method\E\s+(.*)/) {
      print "# id: $query_id, tree algorithm: $method \n" if($verbose);
      print "$1 \n";
      exit;
    }
  }

  die("Found no  $method  tree with id $query_id in it.\n");
}

if (defined $n) {
  my $qid = '';
  while (<>) {
    $qid = $1 if(/^Id\s+(\S+)/);
    if (/^\Q$method\E\s+(.*)/ ) {
      $n_trees_read++;
      if ($n_trees_read == $n) {
	print "# $n" . "th tree using tree algorithm: $method. id: $qid \n" if($verbose);
	print "$1 \n";
	exit;
      }
    }
  }
  die($n . "th tree requested, but only $n_trees_read trees in file.\n");
}
