#!/usr/bin/perl -w
use strict;

my %ids = ();

my $idfile = shift;

open my $fh, "<", "$idfile" or die "couldnt open $idfile for reading.\n";
while (<$fh>) {
  next if (/^\s*#/);
  my @cols = split(" ", $_);
  my $id = shift @cols;
  $id .= '.1';
  #  print "X: [$id] \n";
  $ids{$id}++;
}


while (my $the_line = <>) {
  if($the_line =~ /^\s*#/){
    print $the_line;
next;
  }
  if ($the_line =~ /^\s*(\S+)/) {
    my $the_id = $1;
    #   print "Y: [$the_id] \n";

    my @cols = split(" ", $the_line);
    my $sp = "   ";
    #   print "XXX: ", join(", ", @cols), "\n";
    #exit;
    my $s = $cols[0] . $sp . $cols[1] . $sp .
      join(" ", @cols[2..5]) . $sp .
	join(" ", @cols[6..9]) . $sp .
	  join(" ", @cols[10..13]) . $sp .
	    join(" ", @cols[14..17]) . $sp .
	      join(" ", @cols[18..21]);
my $nest1in3 = ($cols[2] < $cols[10])? 1:0;
my $nest3in5 = (($cols[18] < 0)  or  ($cols[10] < $cols[18]))? 1:0;
my $ndisin3 = $cols[12];
my $OK = ($nest1in3 and $nest3in5 and ($ndisin3 == 0))? '1': '0';
$nest1in3 = ($nest1in3)? '1':'0';
$nest3in5 = ($nest3in5)? '1':'0';

my $xx = "   $nest1in3 $nest3in5 $ndisin3   $OK ";
    print "$s $xx \n" if(exists $ids{$the_id});
  }
}
