#!/usr/bin/perl -w
use strict;

# e.g.:
# KS_test.pl fam9877.nex  1  2000
# read in two files *.run1.p and *.run2.p
# find Kolmogorov-Smirnov 'D' for data in
# column $col (0-based, default = 1)
# skip pts with ngen < 2000

my $bigneg = -1e300;

my $file_basename = shift;
my $datacol = shift; # data column to use.
$datacol = 1 if(!defined $datacol);
my $ngen_skip = shift;
$ngen_skip = 3000 if (!defined $ngen_skip);

# store data in hashes
my @val_count_hashes = ({}, {}); #
my @counts = (0, 0);
foreach my $i (0..1){
	my $irun = $i+1;
	my $filename = "$file_basename.run" . $irun . ".p"; 
#print "filename: $filename\n";
	open my $fh, "<$filename";
	while(<$fh>){
		my @cols = split(" ", $_);
# skip non-numerical stuff.
		next unless($cols[0] =~ /^\d+$/); 
		my ($ngens, $x) = @cols[0,$datacol];
#	print "ngens x: $ngens  $x\n";		
		next if($ngens < $ngen_skip);
		$val_count_hashes[$i]->{$x}++;
		$counts[$i]++;
	}
#	print "$i $counts[$i] \n";
	close $fh;
}

# get cumulative distributions:
my @val_cumeprob_hashes = ({$bigneg=>0}, {$bigneg=>0});
foreach my $i (0..1){
	my $cume_prob = 0;
	foreach (sort {$a <=> $b} keys %{$val_count_hashes[$i]}){
		$cume_prob += $val_count_hashes[$i]->{$_}/$counts[$i];
		$val_cumeprob_hashes[$i]->{$_} = $cume_prob;
#		print "$i $_ $cume_prob ", $val_count_hashes[$i]->{$_}, "\n";	
	}
#	print "\n";
}

print "K-S D: ", Kolmogorov_Smirnov_D(@val_cumeprob_hashes), "\n";

sub Kolmogorov_Smirnov_D{
# get the maximum difference between two empirical cumulative distributions.
# Arguments are two hashrefs, each representing an empirical cumulative distribution.
# Each key is a data point (a real number), and the corresponding hash value is
# the proportion of data pts <= to it. So the largest key should have value 1.
	my $val_cumeprob1 = shift; 
	my $val_cumeprob2 = shift;
	my @sorted_vals1 = sort {$a <=> $b} keys %{$val_cumeprob1};
	my @sorted_vals2 = sort {$a <=> $b} keys %{$val_cumeprob2};
#	print "sorted vals1: ", join(",", @sorted_vals1), "\n\n";
#  print "sorted vals2: ", join(",", @sorted_vals2), "\n\n";

	my ($i1, $i2) = (0, 0);

	my $size1 = scalar @sorted_vals1;
	my $size2 = scalar @sorted_vals2;

	my $D = 0; 
	while(1){
#	print "$i1 $i2 \n";
		my ($xlo1, $xhi1) = ($sorted_vals1[$i1], $sorted_vals1[$i1+1]);
		my ($xlo2, $xhi2) = ($sorted_vals2[$i2], $sorted_vals2[$i2+1]);
#			print "$xlo1, $xhi1   $xlo2, $xhi2 \n";
		die "$xlo1 > $xhi2 ??\n" if($xlo1 > $xhi2);
		die "$xlo2 > $xhi1 ??\n" if($xlo2 > $xhi1);

		my ($cume_prob1, $cume_prob2) = ($val_cumeprob1->{$xlo1}, $val_cumeprob2->{$xlo2});
		my $abs_diff = abs($cume_prob1 - $cume_prob2);
		$D = $abs_diff if($abs_diff > $D);
		if($xhi1 <= $xhi2){
			$i1++;
		}elsif($xhi2 <= $xhi1){
			$i2++;
		}else{
			die "$xhi1 xhi2 should be numerical.\n";
		}
		last if($i1 == $size1-1);
		last if($i2 == $size2-1);
	}
	return $D;
}




