#!/usr/bin/perl -w
use strict;

# filter a  *.fastas (or *.alfastas) file,
# keeping only the families which have fasta sequences present.

my $print_this = 0;
while(my $the_line = <>){
	if($the_line =~ /^Id (\S+)/){
		my $next_line = <>;
		if($next_line =~ /^\s+$/){ # blank line - don't print this seq.
			$print_this = 0;
		}elsif($next_line =~ /^>\S+/){ # fasta present - start printing (including the 'Id Medtr...' line) 
			$print_this = 1;
			print $the_line;
			print $next_line;
		}else{
			warn "Line after Id Medtr... line is of unexpected format: \n" . $the_line;
		}
	}else{
		if($print_this){
			print $the_line;
		}
	}
}
