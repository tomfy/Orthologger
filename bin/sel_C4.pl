#!/usr/bin/perl
use strict;

while (<>) {
  my @cols = split(" ", $_);
  my $id = shift @cols;
  my $type = shift @cols;
  # heights above query of various subtree root nodes:
  my ($C4M3of4, $C4M4of4, $C3M, $C3D, $Basal) = ( $cols[0], $cols[5], $cols[10], $cols[15], $cols[20]); 

  my $C4M = $C4M4of4;


  #  print "$C4M  $C3M  $C3D  $Basal \n";
  #     print;

  next if($C4M < 0);
  if ($Basal < 0) { # no basals
    if ($C3M < 0  and  $C3D < 0) { # no C3s. Tree has just the 4 C4's
 #     print "$C4M  $C3M  $C3D  $Basal \n";
      print;
      #exit;
    }
  } else {
    if ($Basal > $C4M) { # basal present. 
      if ( ($C3M < 0  or  $C3M > $Basal)  and  ($C3D < 0 or $C3D > $Basal) ) { # must go higher than first basal to get (other) C3s
#	print "$C4M  $C3M  $C3D  $Basal \n";
	print;
      }
    }
  }
}

