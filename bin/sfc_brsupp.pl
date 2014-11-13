#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);

my $min_support = shift || 0.0;

my $require_no_amborella_in_monocots = shift || 0;

my ($AMpos_offset, $monocots_offset, $basals_offset, $non_angiosperms_offset, $AMneg_offset) = 
  (0, 5, 10, 15, 20);

while(<>){
	next if(/^\s*#/);
	next if(/^\s*$/);

	my @cols = split(" ", $_);
	my $id = shift @cols;
	my $type = shift @cols;
	if(/^Medtr/){
#	  print STDERR  join(", ", @cols), "\n";
		next if($cols[$AMpos_offset] <= 0); # no AMpos dicots clade
#print STDERR "AAAA\n";
		next if($cols[$monocots_offset] <= 0); # no monocots clade
#print STDERR "BBBB\n";
		next if( ($cols[$monocots_offset] > 0) and ($cols[$monocots_offset] <= $cols[$AMpos_offset]) ); # 6/12D not nested in 3/6M
#print STDERR "CCCC\n";
		if($require_no_amborella_in_monocots){
	       	next if( ($cols[$basals_offset] > 0) and ($cols[$basals_offset] <= $cols[$monocots_offset]) ); # basals present in 3/6M
		}
	next if( ($cols[$non_angiosperms_offset] > 0) and ($cols[$non_angiosperms_offset] <= $cols[$monocots_offset]) ); # basals present in 3/6M
# print STDERR "XXX: ", $cols[5], "  ", $cols[15], "\n";
			next if( ($cols[$AMneg_offset] > 0) and ($cols[$AMneg_offset] <= $cols[$monocots_offset]) ); # negatives present in 3/6M
# otherwise - it's OK
#	  print "EEEE $_ \n";
			my $Pdicots_in_Monocots_support = max( (1 - ($cols[$monocots_offset+1]/$cols[$AMpos_offset+1])), -1);
		my $Monocots_in_Basals_support = max( (($cols[$basals_offset] > 0)? (1 - ($cols[$basals_offset+1]/$cols[$monocots_offset+1])) : 1), -1);
	my $Monocots_in_Nonangiosperms_support = max( (($cols[$non_angiosperms_offset] > 0)? (1 - ($cols[$non_angiosperms_offset+1]/$cols[$monocots_offset+1])) : 1), -1);
		my $Negatives_not_in_Monocots_support = max( (($cols[$AMneg_offset] > 0)? (1 - ($cols[$AMneg_offset+1]/$cols[$monocots_offset+1])) : 1), -1);
		chomp;
	my $min_supp1 = min($Pdicots_in_Monocots_support, $Monocots_in_Basals_support, $Negatives_not_in_Monocots_support);
		$min_supp1 = max($min_supp1, -1.0);
		my $min_supp2 = min($Pdicots_in_Monocots_support, $Monocots_in_Nonangiosperms_support, $Negatives_not_in_Monocots_support);
#		print STDERR "min_supp: $min_supp1 $min_supp2    $min_support.\n";
		if($min_supp2 >= $min_support){
			print "$_  ";
			printf("%6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f\n", $Pdicots_in_Monocots_support, $Monocots_in_Basals_support, $Monocots_in_Nonangiosperms_support, $Negatives_not_in_Monocots_support, $min_supp1, $min_supp2);
		}
	}
}

# the end
