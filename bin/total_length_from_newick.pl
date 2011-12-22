#!/usr/bin/perl -w

my $newick = <>; # newick string with branch lengths.
 
my $total = 0;
while($newick =~ s/:(-?\d+[.]\d+)//){

my $branch_length = $1;
$total += $branch_length;
print "bl: $branch_length  total: $total \n";
}
print $total, "\n";
