#!/usr/bin/perl -w
use strict;

# m8 -> abc
# impose max e-value (i.e. keep only matches with e-value <= $max_eval)
# keep only best solution for each id1-id2 pair


my $max_eval = shift || 0.00000001; # 1e-8 default.

my $old_idpair = 'xxxxxxx xxxxxxx' ;
my $the_line = '';
while($the_line = <>){
	my @cols = split(" ", $the_line);
	my ($id1, $id2, $evalue) = @cols[0,1,10];
	my $idpair = "$id1 $id2";

	if(($idpair ne $old_idpair)  and  ($evalue <= $max_eval)){
		print "$idpair $evalue \n";  #$the_line;
		$old_idpair = $idpair;
	}
}
