package TlyUtil;
use strict;
use List::Util qw ( min max sum );

# 


sub Kolmogorov_Smirnov_D{
# get the maximum difference between two empirical cumulative distributions.
# Arguments are two hashrefs, each representing an empirical cumulative distribution.
# Each key is a data point (a real number), and the corresponding hash value is
# the proportion of data pts <= to it. So the largest key should have value 1.
	my $val_cumeprob1 = shift; 
	my $val_cumeprob2 = shift;
	my @sorted_vals1 = sort {$a <=> $b} keys %{$val_cumeprob1};
	my @sorted_vals2 = sort {$a <=> $b} keys %{$val_cumeprob2};

	my ($i1, $i2) = (0, 0);

	my $size1 = scalar @sorted_vals1;
	my $size2 = scalar @sorted_vals2;

	my $D = 0;
	while(1){
		my ($xlo1, $xhi1) = ($sorted_vals1[$i1], $sorted_vals1[$i1+1]);
		my ($xlo2, $xhi2) = ($sorted_vals2[$i2], $sorted_vals2[$i2+1]);
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

1;


