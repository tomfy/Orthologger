#!/usr/bin/perl -w
use strict;

my $id1;
my $id2;
my $e_value;

while(<>){
   next if(/^#/);
# if(/^\[\s*(\S+)\s*\]/){ # older format: [ id1 id2 id3 ], (new is id1:id2:id3 )
if(/^\s*(\S+)/){
   $id1 = $1;
}elsif(/^\s+(\S+)\s+(\S+)/){
   $id2 = $1;
   $e_value = $2;
   print "$id1  $id2  $e_value \n";
}elsif(/^\s*$/){
# do nothing - blank line in input
}else{
   print STDERR "other kind of line: {$_} \n";
   exit;
}
}

