#!/usr/bin/perl -w
use strict;


my $the_col = shift // 10;

my $max_eval_so_far = -1;
my %id_eval = ();
while (<>) {
 # print;
  next if(/^\s*$/); next if(/^\s*#/);
  my @cols = split(" ", $_);
  my $evcol;
  if (scalar @cols == 12) {	# m8 format
    $evcol = 10;
  } elsif (scalar @cols == 3) { # abc format
    $evcol = 2;
  } else {
    die "Number of cols is ", scalar @cols , " format is ???\n";
  }
#print $evcol, "  ", join(";  ", @cols), "\n";
  my ($id1, $id2, $eval) = @cols[0,1,$evcol];
#print "$id1 $id2 $eval $max_eval_so_far \n";
#exit;
  $id_eval{$id1} = $eval;
  if ($eval > $max_eval_so_far) {
    $max_eval_so_far = $eval;
   print "Max eval so far: $max_eval_so_far   $id1 $id2\n";
  }
}
exit;
my $log10 = log(10);
for (keys %id_eval) {
  my $ev = $id_eval{$_};
  print "$_   $ev   ", ($ev > 0)? log($ev)/$log10 : -10000, "\n";
}
print "Max eval: $max_eval_so_far \n";
