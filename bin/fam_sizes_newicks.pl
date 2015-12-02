#!/usr/bin/perl -w
use strict;

my ($id, $fs, $numb_eq);

while(<>){ if(/^FT /){ 
   $numb_eq = $_ =~ tr/=/=/; 
   #print $numb_eq, "\n";
}elsif(/^Id (\S+)\s.*fam_size:\s+(\S+)/){
   $id = $1; $fs = $2;
}elsif(/^\s*$/){
   print "$id  $fs  $numb_eq \n" if($fs != $numb_eq);
}
}
