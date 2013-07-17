#!/usr/bin/perl -w
use strict;

my $idpattern = shift;
my $print = 0;
while(<>){
	if(/^Id /){
		if(/$idpattern/){
			$print = 1;
#			print "$_ turning on print \n";
		}else{
			last if($print == 1); 
		}	
#		$print = 0;
	}else{
		print $_ if($print);
	}
}

