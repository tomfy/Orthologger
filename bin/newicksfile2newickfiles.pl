#!/usr/bin/perl -w
use strict;

my $newicks_file = shift;

open my $fhin, "<", "$newicks_file" or die "Couldn't open $newicks_file for reading.\n";

while(<$fhin>){
if(/^Id\s+(\S+)/){
my $out_filename = $1 . ".newick";
open my $fhout, ">", $out_filename;
my $newick_string = <$fhin>;
$newick_string =~ s/^[^(]+//; # delete everything before first '('
print $fhout $newick_string;
close $fhout;
}

}
