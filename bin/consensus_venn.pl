#!/usr/bin/perl -w
use strict;

# take a pair cladesout files as input.
# define some criteria for categories, such as
# 'D<M<S, no M paralogs, no Brassicas in S'
# and require consensus of both input files for an id to
# be put in the category

my $in_file_1 = shift;		# e.g. based on muscle alignments
my $in_file_2 = shift;		# e.g. based on mafft alignments

my @select_option_strings = (
			     " -nest '1<3' -maxdisallowed '3:0' ", # outer: D(7/8) < M(3/4), No Brassicas in M

			     " -nest '1<3,3<4' -maxdisallowed '4:0' ", # D<M<S, No Brassicas in S
			     " -nest '1<3' -maxpara '3:0' -maxdisallowed '3:0' ", # D<M, , No paralogs in M, No Brassicas in M
			     " -nest '1<3,3<4' -maxpara '3:0' -maxdisallowed '4:0' " # innermost region: D<M<S, no paralogs in M, no Brassicas in S
			    );

 @select_option_strings = (
 			     " -nest '2<3' -maxdisallowed '3:0' ", # outer: D(7/8) < M(3/4), No Brassicas in M

 			     " -nest '2<3,3<4' -maxdisallowed '4:0' ", # D<M<S, No Brassicas in S
 			     " -nest '2<3' -maxpara '3:0' -maxdisallowed '3:0' ", # D<M, , No paralogs in M, No Brassicas in M
 			     " -nest '2<3,3<4' -maxpara '3:0' -maxdisallowed '4:0' " # innermost region: D<M<S, no paralogs in M, no Brassicas in S
 			    );
my %id_either = (); # ids which are OK in at least one of the input cladesout files
my @Ns = ();
my @id_hashes = ();
for my $sel_opt_str (@select_option_strings) {
  my %cons_id_hash = ();
  my ($N1, $N2) = (0, 0);
  my $count_acc = 0;
  my %id_accinfo_1 = ();
  # print "$sel_opt_str \n";
  my $cl1 = " select_from_cladesout.pl $sel_opt_str < $in_file_1 "; # | grep -P '^\s*(\S+\s+){2}1' ";
  for ( split("\n", `$cl1`) ) {
    /^\s*(\S+)\s+(.*)/;
    my $acc_info_1 = [split(" ", $2)];
    $id_accinfo_1{$1} = $acc_info_1;
    if ($acc_info_1->[1] == 1) {
      $N1++;
      $id_either{$1} += 10000;
    }
    # print "[", join("; ", @{$id_accinfo_1{$1}}), "]\n";
  }
  my $cl2 = " select_from_cladesout.pl $sel_opt_str < $in_file_2 "; # | grep -P '^\s*(\S+\s+){2}1' ";
  for ( split("\n", `$cl2`) ) {
    /^\s*(\S+)\s+(.*)/;
    my $id2 = $1;
    my $acc_info_1 = (exists $id_accinfo_1{$id2})?  $id_accinfo_1{$id2}: [0, 0, 0, 0];
    my $acc_info_2 = [split(" ", $2)]; # $2;
    #print "$id2  ", $acc_info_1->[1], "  ", $acc_info_2->[1], "\n";
    if ($acc_info_2->[1] == 1) {
      $N2++;
      $id_either{$id2}++;
      if ($acc_info_1->[1] == 1) {
	#     print "$id2   ", join(",", @$acc_info_1), "    ", join(",", @$acc_info_2), "\n";
	$count_acc++;
	$cons_id_hash{$id2}++;
      }
    }
  }
  push @Ns, $count_acc;
  print "# Nacc: $count_acc  $sel_opt_str \n";
  print "# N1, N2, Nunion, Nintersection: $N1  $N2 ", $N1 + $N2 - $count_acc, "  $count_acc \n";
  push @id_hashes, \%cons_id_hash;
}

my $Nouter =  $Ns[0] - ($Ns[1] + $Ns[2]) + $Ns[3];
my $Nupper = $Ns[1] - $Ns[3];
my $Nmiddle = $Ns[3];
my $Nlower = $Ns[2] - $Ns[3];

print "# Outer, upper, middle, lower: $Nouter  $Nupper  $Nmiddle  $Nlower \n";

venn($id_hashes[0], $id_hashes[1], $id_hashes[2]  , \%id_either);

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
