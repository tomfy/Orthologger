#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );

while (<>) {
  chomp;
  s/;$//;
  my $ordered_newick_string = order_newick($_);
  print $ordered_newick_string, "\n";
}


sub order_newick{
  my $newick = shift;
# print "top of order_newick. newick: $newick\n";
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
      if ($_ eq '('){ $lmr_paren_count++; }
      if ($_ eq ')'){ $lmr_paren_count--; }
#      print "il,ir, curr char, paren count: $il $ir $_ $n_chars $lmr_paren_count [",
#	join('',@newick_chars[$il..$ir-1]), "]\n";
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
    }  # loop over chars in @newick_chars
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
