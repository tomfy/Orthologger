#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# take a pair cladesout files as input.
# define some criteria for categories, such as
# 'D<M<S, no M paralogs, no Brassicas in S'
# and require consensus of both input files for an id to
# be put in the category



# Note: in the following array, elements 1 and 2 should specify subsets of elem. 0, and elem. 3 should specify 
# the intersection of 1 and 2. 
my @select_option_strings = ( # 7 of 8 dicots
			     " -nest '1<3' -maxdisallowed '3:0' ", # outer: D(7/8) < M(3/4), No Brassicas in M
			     " -nest '1<3,3<4' -maxdisallowed '4:0' ", # D<M<S, No Brassicas in S
			     " -nest '1<3' -maxpara '3:0' -maxdisallowed '3:0' ", # D<M, , No paralogs in M, No Brassicas in M
			     " -nest '1<3,3<4' -maxpara '3:0' -maxdisallowed '4:0' " # innermost region: D<M<S, no paralogs in M, no Brassicas in S
			    );

 #      @select_option_strings = ( # 6 of 7 dicots
 #      			     " -nest '2<3' -maxdisallowed '3:0' ", # outer: D(6/7) < M(3/4), No Brassicas in M
 #      			     " -nest '2<3,3<4' -maxdisallowed '4:0' ", # D<M<S, No Brassicas in S
 #      			     " -nest '2<3' -maxpara '3:0' -maxdisallowed '3:0' ", # D<M, , No paralogs in M, No Brassicas in M
 #      			     " -nest '2<3,3<4' -maxpara '3:0' -maxdisallowed '4:0' " # innermost region: D<M<S, no paralogs in M, no Brassicas in S
 #      			    );

 # @select_option_strings = ( # 7 of 8 dicots D7/8 < M1/4
 # 			     " -nest '1<3' -maxdisallowed '5:0' ", # outer: D(7/8) < M(3/4), No Brassicas in M(3<4)
 # 			     " -nest '1<3,5<6' -maxdisallowed '6:0' ", # D<M<S, No Brassicas in S
 # 			     " -nest '1<3' -maxpara '5:0' -maxdisallowed '5:0' ", # D<M, , No paralogs in M, No Brassicas in M
 # 			     " -nest '1<3,5<6' -maxpara '5:0' -maxdisallowed '6:0' " # innermost region: D<M<S, no paralogs in M, no Brassicas in S
 # 			    );

 # @select_option_strings = ( # 7 of 8 dicots D7/8 < M2/4
 #  			     " -nest '1<4' -maxdisallowed '5:0' ", # outer: D(7/8) < M(3/4), No Brassicas in M(3<4)
 #  			     " -nest '1<4,5<6' -maxdisallowed '6:0' ", # D<M<S, No Brassicas in S
 #  			     " -nest '1<4' -maxpara '5:0' -maxdisallowed '5:0' ", # D<M, , No paralogs in M, No Brassicas in M
 #  			     " -nest '1<4,5<6' -maxpara '5:0' -maxdisallowed '6:0' " # innermost region: D<M<S, no paralogs in M, no Brassicas in S
 #  			    );

  # @select_option_strings = ( # 7 of 8 dicots D7/8 < M3/4
  # 			     " -nest '1<5' -maxdisallowed '5:0' ", # outer: D(7/8) < M(3/4), No Brassicas in M(3<4)
  # 			     " -nest '1<5,5<6' -maxdisallowed '6:0' ", # D<M<S, No Brassicas in S
  # 			     " -nest '1<5' -maxpara '5:0' -maxdisallowed '5:0' ", # D<M, , No paralogs in M, No Brassicas in M
  # 			     " -nest '1<5,5<6' -maxpara '5:0' -maxdisallowed '6:0' " # innermost region: D<M<S, no paralogs in M, no Brassicas in S
  # 			    );


 
my $in_file_1 = undef;		# e.g. based on muscle alignments
my $in_file_2 = undef;		# e.g. based on mafft alignments
my $min_sum_acc_bs = 0; # default minimum number of accepted bootstraps (sum of analyses 1 and 2)
my $min_each_acc_bs = 0; # default min number of accepted bootstraps (must be >= this in both analysis 1 and 2)

GetOptions(
    'file1=s' => \$in_file_1,
    'file2=s' => \$in_file_2,
    'min_sum_acc_bs=i' => \$min_sum_acc_bs,
    'min_both_acc_bs=i' => \$min_each_acc_bs,
);
die "input file 1, ($in_file_1), undefined or no such file.\n" if(!defined $in_file_1  or ! -f $in_file_1);
die "input file 2, ($in_file_2), undefined or no such file.\n" if(!defined $in_file_2  or ! -f $in_file_2);
print "# Input files: $in_file_1, $in_file_2 \n";
print "# both n_accepted bootstraps >= $min_each_acc_bs; sum n_accepted bootstraps >= $min_sum_acc_bs \n";
my %id_either = (); # ids which are OK in at least one of the input cladesout files
my @Ns = ();
my @id_hashes = ();
for(my $ii=0; $ii < scalar @select_option_strings; $ii++){
 
my $sel_opt_str = $select_option_strings[$ii];
# print $sel_opt_str, "\n";
  my %cons_id_hash = ();
  my ($N1, $N2, $N1and2) = (0, 0, 0);
  my $count_acc = 0;
  my %id_accinfo_1 = ();
  # print "$sel_opt_str \n";
  my $cl1 = " select_from_cladesout.pl $sel_opt_str < $in_file_1 "; # | grep -P '^\s*(\S+\s+){2}1' ";
#print "$cl1 \n";
  for ( split("\n", `$cl1`) ) {
#    print "sfc output line: $_ \n";
    next if(/^\s*#/  or  /^\s*$/);
    /^\s*(\S+)\s+(.*)/;
#print "1: $1    2: $2 \n";
    my $acc_info_1 = [split(" ", $2)];
#    print "acc_info_1:  ", join(";", @$acc_info_1), " \n";
    $id_accinfo_1{$1} = $acc_info_1;
    if ($acc_info_1->[1] == 1) {
      $N1++;
      $id_either{$1} += 10000;
    }
    # print "[", join("; ", @{$id_accinfo_1{$1}}), "]\n";
  }
  my $cl2 = " select_from_cladesout.pl $sel_opt_str < $in_file_2 "; # | grep -P '^\s*(\S+\s+){2}1' ";
  for ( split("\n", `$cl2`) ) {
#    print "$_\n";
  next if(/^\s*#/  or  /^\s*$/);
#    print "$_ \n";
    /^\s*(\S+)\s+(.*)/;
    my $id2 = $1;
#  print "file2, id: $1,  info: $2 \n";
    my $acc_info_1 = (exists $id_accinfo_1{$id2})?  $id_accinfo_1{$id2}: [0, 0, 0, 0];
    my $acc_info_2 = [split(" ", $2)]; # $2;
#    print "$id2  ", $acc_info_1->[1], "  ", $acc_info_2->[1], "\n"; exit;
    if ($acc_info_2->[1] == 1) {
      $N2++;
      $id_either{$id2}++;
      if ($acc_info_1->[1] == 1) {
	$N1and2++;
	#     print "$id2   ", join(",", @$acc_info_1), "    ", join(",", @$acc_info_2), "\n";
	if(($acc_info_1->[3] >= $min_each_acc_bs and $acc_info_2->[3] >= $min_each_acc_bs) and ($acc_info_1->[3] + $acc_info_2->[3] >= $min_sum_acc_bs)){
	$count_acc++;
	$cons_id_hash{$id2}++;
      }
      }
    }
  }
  push @Ns, $count_acc;
  print "# Region $ii;   ", $select_option_strings[$ii], "\n";
  print "# Nfile1, Nfile2, Nunion, Nintersection: $N1  $N2 ", $N1 + $N2 - $N1and2, "  $N1and2 (No requirement on bootstraps)\n";
  print "# Nintersection and satisfying bootstrap requirements: $count_acc  $sel_opt_str  $min_each_acc_bs  $min_sum_acc_bs \n\n";
  push @id_hashes, \%cons_id_hash;

}

my $Nouter =  $Ns[0] - ($Ns[1] + $Ns[2]) + $Ns[3];
my $Nupper = $Ns[1] - $Ns[3];
my $Ninner = $Ns[3];
my $Nlower = $Ns[2] - $Ns[3];
my $Check = $Nupper + $Nlower + 2*$Ninner - ($Ns[1] + $Ns[2]);


my $cons_id_hashref = $id_hashes[0]; # consensus least restrictive region
for(keys %$cons_id_hashref){
  my ($in1, $in2, $in3) = ($id_hashes[1]->{$_}, $id_hashes[2]->{$_}, $id_hashes[3]->{$_});
  my $region = '1 0 0 0 outer';
  if($in3){
    $region = '1 1 1 1 inner';
  }elsif($in1){
    $region = '1 1 0 0 upper crescent';
  }elsif($in2){
    $region = '1 0 1 0 lower crescent';
  }

  print "$_   $region\n";
}

print "\n";
print "# both n_accepted bootstraps >= $min_each_acc_bs; sum n_accepted bootstraps >= $min_sum_acc_bs \n";
print "# 1 0 0 0 Outer:          $Nouter    (In 0 but not in 1, 2, or 3) \n" . 
  "# 1 1 0 0 Upper crescent: $Nupper    (In 1 but not 2) \n" . 
  "# 1 1 1 1 Inner:          $Ninner    (In both 1 and 2) \n" . 
  "# 1 0 1 0 Lower crescent: $Nlower    (In 2 but not 1) \n";

# , upper, middle, lower: $Nouter  $Nupper  $Ninner  $Nlower \n";
exit;
venn($id_hashes[0], $id_hashes[1], $id_hashes[2]  , \%id_either);
#exit;
for(keys %id_either){
  print "$_    notboth  ", $id_either{$_}, "\n";
}



sub venn{
my $id_hash1 = shift;
my $id_hash2 = shift;
my $id_hash3 = shift;
my $id_either = shift;
#my $id_hash4 = shift;

for my $id1 (keys %$id_hash1){
  delete $id_either->{$id1};
  if(exists $id_hash2->{$id1}){
   if(exists $id_hash3->{$id1}){
      print "$id1  middle \n";
    }else{
      print "$id1  upper \n";
    }
  }else{
    if(exists $id_hash3->{$id1}){
      print "$id1  lower \n";
    }else{
      print "$id1  outer \n";
    }
  }
}

}

sub union{
  my $id_hash1 = shift;
  my $id_hash2 = shift;
my @union_ids = ();
  for my $id1 (keys @$id_hash1){
    if(exists $id_hash2->{$id1}){
      push @union_ids, $id1;
}
  }
  return \@union_ids;
}
