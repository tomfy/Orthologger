package TlyUtil;
use strict;
use List::Util qw ( min max sum );

my $bigneg = -1e300;

# given a set of numbers (some of which may occur more than once), which
# are stored as keys in a hash, with the values being how many times they occur
# or more generally whatever weights you want, sort the keys and get the
# cumulative distribution.
sub cumulative_prob{
  my $val_weight_href = shift; # hashref, key: numbers, values: weights.
  my $sum_weights = shift;
  my $val_cumeprob_href = { $bigneg => 0};
  my $cume_prob = 0;
  foreach (sort {$a <=> $b} keys %{$val_weight_href}) {
    $cume_prob += $val_weight_href->{$_}/$sum_weights;
    $val_cumeprob_href->{$_} = $cume_prob;
    #		print "$i $_ $cume_prob ", $val_count_hashes[$i]->{$_}, "\n";	
  }
  #	print "\n";
  return $val_cumeprob_href;
}


sub Kolmogorov_Smirnov_D{
# get the maximum difference between two empirical cumulative distributions.
# Arguments are two hashrefs, each representing an empirical cumulative distribution.
# Each key is a data point (a real number), and the corresponding hash value is
# the proportion of data pts <= to it. So the largest value should be 1.
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

# operates on a newick of form (3,(6,4))
# i.e. no whitespace, no branch lengths, ids must be numbers.
# so just parens, commas and numbers
# puts the leaves in order, such that at each node the
# subtree with smaller value is on left. The value of an
# internal node is the min of the values of the two child
# nodes, and the value of a leave is its id, which must be a number.
sub order_newick{
  my $newick = shift;
  if ($newick =~ /^(\d+)$/) {	# subtree is leaf!
    return ($1, $newick);
  } else {			# subtree has > 1 leaf.
    my %label_newick = ();
    $newick =~ /^[(](.*)[)]$/;
    my @newick_chars = split('',$1); # without surrounding ()
    my $lmr_paren_count = 0;
    my ($il, $ir) = (0, 0);
    my $n_chars = scalar @newick_chars;
    my $min_label = 10000000;
    foreach (@newick_chars) {
      die "$_ ", $newick_chars[$ir], " not same!\n" if($_ ne $newick_chars[$ir]);
      if ($_ eq '(') {
	$lmr_paren_count++;
      }
      if ($_ eq ')') {
	$lmr_paren_count--;
      }

      if (($ir == $n_chars-1) or ($_ eq ',' and $lmr_paren_count == 0)) { #split
	my $ilast = ($ir == $n_chars-1)? $ir : $ir-1;
	my $sub_newick = join('', @newick_chars[$il..$ilast]);
	#       print "subnewick $sub_newick\n";
	my ($label, $ordered_subnewick) = order_newick($sub_newick);
	$label_newick{$label} = $ordered_subnewick;
	$min_label = min($min_label, $label);
	$il = $ir+1; $ir = $il; # skip the ','
      } else {
	$ir++;
      }
    }				# loop over chars in @newick_chars
    my $ordered_newick = '';
    foreach (sort {$a <=> $b} keys %label_newick) {
      $ordered_newick .= $label_newick{$_} . ",";
    }
    $ordered_newick =~ s/,$//;
    $ordered_newick = '(' . $ordered_newick . ')';
    return ($min_label, $ordered_newick);
  }
  die "shouldnt get here, in order_newick\n";
}


1;


