#!/usr/bin/perl -w
use strict;

my $old_idpair = 'xxxxxxx xxxxxxx' ;
my $the_line = '';
while($the_line = <>){
	my @cols = split(" ", $the_line);
	my ($id1, $id2) = @cols[0,1];
	my $idpair = "$id1 $id2";

	if($idpair ne $old_idpair){
		print $the_line;
		$old_idpair = $idpair;
	}
}
