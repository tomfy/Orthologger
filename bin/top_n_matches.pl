#!/usr/bin/perl -w
use strict;

# keep the top $max_matches matches for each query (m8 or abc formats)

my $max_matches = shift || 100; 

my $old_id1 = 'xxxxxxx_xxxxxxx' ;
my $the_line = '';
my $count = 0;
while($the_line = <>){
	my @cols = split(" ", $the_line);
	my ($id1, $id2) = @cols[0,1];

	if($id1 ne $old_id1){
		$old_id1 = $id1;
		$count = 1;
	}else{
		$count++;
	}
#	print "$id1  $count \n";
	if($count <= $max_matches){
		print $the_line;
	}
}
