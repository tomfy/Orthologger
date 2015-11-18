#!/usr/bin/perl -w
use strict;

my $fam_count = 0;
my $fam_size = -100;
my ($idline_famsize, $old_idline_famsize) = (0, 0);
my ($old_id, $id) = ('no_id', 'no_id');
while(<>){
	if(/^Id\s+(\S+).*fam_size:\s+(\d+)/){
   $old_idline_famsize = $idline_famsize;
   $old_id = $id;
   $id = $1;
   $idline_famsize = $2;
	print "$old_id   $fam_count   $fam_size   $old_idline_famsize\n" if($fam_size != $old_idline_famsize);
	$fam_size = 0;
	$fam_count++;
}elsif(/^>/){
	$fam_size++;
}
}
print "$id   $fam_count   $fam_size   $idline_famsize\n" if($fam_size != $idline_famsize);
print "# $fam_count families. Done.\n";
