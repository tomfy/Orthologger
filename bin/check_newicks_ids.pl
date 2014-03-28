#!/usr/bin/perl -w
use strict;

# to check a 'newicks' format file (i.e. output of fttree.pl)
# to see that query ids are present in the corresponding trees

my ($q_count, $tree_count, $OK_count, $bad_count) = (0,0,0,0);

while (<>) {
    if (/^Id (Med\S+)/) {
      $q_count++;
        my $id     = $1;
        my $newick = <>;
        if ( $newick =~ /^\(/ ) {
	  $tree_count++;
           $newick =~ s/;?\s*$//;
            if ( $newick =~ /\Q$id\E/ ) {
           #     print "$id  OK\n";
		$OK_count++;
            }
            else {
                print "# $id  Not found\n";
		$bad_count++;
            }
        }
    }
}
print "$q_count, $tree_count, $OK_count, $bad_count \n";
