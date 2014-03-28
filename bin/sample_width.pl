#!/usr/bin/perl -w
use strict;

my $dimensions = shift || 3;
my $N = 50000;

for (my $sample_size = 20; $sample_size < $N; $sample_size = int($sample_size*2)) {
  my $count = 0;
  my ($sum1, $sum2, $sum3, $sum10, $sum20) = (0, 0, 0, 0, 0);
  while ($count*$sample_size < $N) {
    my @xs = ();
    for (1..$sample_size) {
      push @xs, rand_gauss(20, $dimensions);
    }
    $sum1 += lnw(\@xs, 1);
    $sum2 += lnw(\@xs, 2);
    $sum3 += lnw(\@xs, 3);
    $sum10 += lnw(\@xs, 10);
    $sum20 += lnw(\@xs, 20);
    $count++;
  }
  print "$sample_size  ", $sum1/$count, "  ", $sum2/$count, "  ", $sum3/$count, "  ", $sum10/$count, "  ", $sum20/$count, "\n";
}

sub rand_gauss{
  my $n = shift || 200;
  my $d = shift || 1;
  my @vector = ();
  for (1..$d) {
    my $sum = 0;
    for (1..$n) {
      $sum += rand() - 0.5;
    }
    push @vector, $sum/$n;
  }
  return \@vector;
}


sub lnw{
  my $xs = shift;
  my $power = shift || 2;
  my $sum = 0.0;
  for my $point (@$xs) {
    my $s = 0;
    for (@$point) {
      $s += $_**2;
    }
    $sum += abs($_)**$power;
  }
  $sum /= scalar @$xs;
  return $sum**(1/$power);
}
