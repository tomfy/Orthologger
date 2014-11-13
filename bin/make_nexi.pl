#!/usr/bin/perl -w
use strict;

my $newicks_filename = shift;

open my $fh, "<", "$newicks_filename" or die "couldn't open [$newicks_filename] for reading.\n";

while(<$fh>){
	if(/^Id (M\S+)\s/){
	my $id = $1;
#print "ID $id \n";
	my $newick_line = <$fh>;
# print $newick_line, "\n";
	$newick_line =~ s/^\s*\S+\s*\(/(/; #remove stuff before first left paren.
#	print $newick_line, "\n";
	open my $fhout, ">", "tmp.newick";
	print $fhout "Id $id \n", "XX ", $newick_line, "\n";
	close $fhout;	
	my $out_nexus_filename = $id . ".nexus";
	system "newicks2figtreenexus.pl -new tmp.newick -id $id > $out_nexus_filename ";
	}
}
