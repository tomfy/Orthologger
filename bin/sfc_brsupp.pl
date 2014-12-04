#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);

my $MinN_min_support = shift || 0.5;
my $DinN_min_support = shift || 0.5;
my $DinM_min_support = shift || 0.5;
my $MinNA_min_support = shift || 0.5;

print "# MinN_min_support, DinN_min_support, DinM_min_support, MinNA_min_support: ",
  "$MinN_min_support, $DinN_min_support, $DinM_min_support, $MinNA_min_support. \n";
my $require_no_amborella_in_monocots = shift || 0;

my ($AMpos_offset, $monocots_offset, $basals_offset, $non_angiosperms_offset, $AMneg_offset) = 
  (0, 5, 10, 15, 20);

while (<>) {
  next if(/^\s*#/);
  next if(/^\s*$/);

  my @cols = split(" ", $_);
  my $id = shift @cols;
  my $type = shift @cols;
  if (/^Medtr/) {
    my $monocots_node_height = $cols[$monocots_offset];
    my $AMposdicots_node_height = $cols[$AMpos_offset];
    my $basals_node_height = $cols[$basals_offset];
    my $nonangiosperms_node_height = $cols[$non_angiosperms_offset];
    my $AMneg_node_height = $cols[$AMneg_offset];

   # next if($AMposdicots_node_height <= 0);	   # no AMpos dicots clade
  #  next if($monocots_node_height <= 0); # no monocots clade

  #  next if($monocots_node_height <= $AMposdicots_node_height); # 6/12D not nested in 3/6M
    # if ($require_no_amborella_in_monocots) {
    #   next if( ($basals_node_height > 0) and ($basals_node_height <= $monocots_node_height) ); # basals present in 3/6M
    # }
    # next if( ($nonangiosperms_node_height > 0) and ($nonangiosperms_node_height <= $monocots_node_height) ); # basals present in 3/6M
    # next if( ($AMneg_node_height > 0) and ($AMneg_node_height <= $monocots_node_height) ); # negatives present in 3/6M

    my $monocots_eps_prod = $cols[$monocots_offset+1];
    my $AMposdicots_eps_prod = $cols[$AMpos_offset+1];
    my $nonangiosperms_eps_prod = $cols[$non_angiosperms_offset+1];
    my $basals_eps_prod = $cols[$basals_offset+1];
    my $AMneg_eps_prod = $cols[$AMneg_offset+1];

    my $Pdicots_in_Monocots_support = a_in_b_support($AMposdicots_eps_prod, $monocots_eps_prod);
    my $Monocots_in_Basals_support =  a_in_b_support($monocots_eps_prod, $basals_eps_prod);
    my $Monocots_in_Nonangiosperms_support =  a_in_b_support($monocots_eps_prod, $nonangiosperms_eps_prod);
    my $Monocots_in_Negatives_support = a_in_b_support($monocots_eps_prod, $AMneg_eps_prod);
    my $Pdicots_in_Negatives_support = a_in_b_support($AMposdicots_eps_prod, $AMneg_eps_prod);

#max( ($AMneg_eps_prod > 0)? 1 - $AMneg_eps_prod/$monocots_eps_prod : 1, -1);

    chomp;

    my $min_supp1 = min($Pdicots_in_Monocots_support, $Monocots_in_Basals_support, $Monocots_in_Negatives_support);
    $min_supp1 = max($min_supp1, -1.0);
    my $min_supp2 = min($Pdicots_in_Monocots_support, $Monocots_in_Nonangiosperms_support, $Monocots_in_Negatives_support);
    #		print STDERR "min_supp: $min_supp1 $min_supp2    $min_support.\n";
 #   if ($min_supp2 >= $min_support) {
      if($Pdicots_in_Monocots_support >= $DinM_min_support and 
	 ($Monocots_in_Negatives_support >= $MinN_min_support and $monocots_node_height != $AMneg_node_height) and
	 $Monocots_in_Nonangiosperms_support >= $MinNA_min_support and
	$Pdicots_in_Negatives_support >= $DinN_min_support){
      print "$_  ";
      printf("%6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f\n", 
	     $Pdicots_in_Monocots_support, $Monocots_in_Basals_support, 
	     $Monocots_in_Nonangiosperms_support, $Monocots_in_Negatives_support, $min_supp1, $min_supp2);
    }
  }
}

sub a_in_b_support{
  my $a_eps_prod = shift;
  my $b_eps_prod = shift;
  if ($a_eps_prod == -1) {
    return 0;
  } else {
    if ($b_eps_prod == -1) {
      return 1;
    } else {
      return a_in_b_support_inner($a_eps_prod, $b_eps_prod);
    }
  }
}

sub a_in_b_support_inner{	    # returns value indicating support for
  # clade a being nested inside of clade b, 
  # based on, e.g. local support values (such as from FastTree) for branches 
  # between roots of a and b
  # 0 < x < 0.5 : b in a preferred,
  # x = 0.5 : indifferent,
  # 0.5 < x <= 1 : a in b preferred.
  my $a_eps_prod = shift;
  my $b_eps_prod = shift;
  my $a_in_b_support;
  if ($b_eps_prod < $a_eps_prod) { # a in b ($b_eps_prod/$a_eps_prod << 1 -> strong support for a in b)
    $a_in_b_support = 1 - 0.5*$b_eps_prod/$a_eps_prod;
  } else { # $a_eps_prod < $b_eps_prod; b in a ($a_eps_prod/$b_eps_prod << 1 -> strong support for b in a)
   
    $a_in_b_support = 0.49999*$a_eps_prod/$b_eps_prod; # 
# print STDERR "XXX $a_eps_prod, $b_eps_prod, $a_in_b_support \n";
  }
  return $a_in_b_support;
}

# the end
