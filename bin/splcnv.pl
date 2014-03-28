#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);

my $first_histogram_line = 1;
my $N = undef;
my $bin_number = 1;
my $sumlogNeff = 0;
my $cslne = 0;
my @lines = <>;
 for my $line (@lines){
	last if($line =~ /^\s*#\s*total/);
	next if ($line =~ /\s*#/);

my @cols = split(" ", $line);
my $label = shift @cols;
my $total = pop @cols;


if($first_histogram_line){
	$N = max(@cols); # Not strictly correct ...
$first_histogram_line = 0;
}
my ($mean, $variance) = mean_variance(\@cols);
my $f_mean = $mean/$N;
my $f_variance = $variance/$N**2;
my $N_eff = ($f_variance > 0)? $f_mean*(1 - $f_mean)/$f_variance : -1;
print $bin_number, "  ", $mean*(1 - $mean/$N)/$N**2, "  ", $variance/$N**2, "  ", $N_eff, "\n";
$bin_number++;
	if($N_eff > 0){
	  $sumlogNeff += log($N_eff);
	  $cslne++;
}
      }
print "N_eff:  ", exp($sumlogNeff/$cslne), "\n";



sub mean_variance{
  my $x = shift;
my ($count, $sum_x, $sum_xsqr) = (scalar @$x, 0, 0);
return (undef, undef) if($count == 0);
  for(@$x){
    $sum_x += $_;
    $sum_xsqr += $_*$_;
  }
my $mean = $sum_x/$count;
my $variance = $sum_xsqr/$count - $mean**2;
return ($mean, $variance);
}



