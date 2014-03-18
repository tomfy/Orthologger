#!/usr/bin/perl -w
use strict;

my $id = shift;
my $print_this_line = 0;
while(<>){
	if(/^Id\s/){  # family id line - turn printing on or off.
		if(/^Id\s+\Q$id\E/){
			$print_this_line = 1;
		}else{
			$print_this_line = 0;
		}
	}else{ # fasta lines - print if it's the right family (i.e. $print_this_line is true)
		print $_ if($print_this_line);
	}
}
