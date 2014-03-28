#!usr/bin/perl -w
use strict;

use lib '/home/tomfy/Orthologger/lib';
use CXGN::Phylo::Cluster;

my $N = shift || 100;

my @labels = (0,1,2,3,4,5);
for(1..$N){
my @xs = map(rand(), @labels);
my @ys = map(rand(), @labels);
my %lp_dist = ();
for (my $i = 0; $i < scalar @labels; $i++) {
  for (my $j = $i+1; $j < scalar @labels; $j++) {
    my $lp = "$i,$j";
    $lp_dist{$lp} = (($xs[$i] - $xs[$j])**2 + ($ys[$i] - $ys[$j])**2)**0.5;
  }
}

my $cluster_obj = CXGN::Phylo::Cluster->new(\@labels, \%lp_dist);

my ($bestness, $string) =  $cluster_obj->best_n_partition(2); #, 'max_inter_min_intra');
print $bestness, "\n";
# for(my $i = 0; $i < scalar @labels; $i++){
#   my $the_label = $labels[$i];
# print "$i  $the_label  ", $xs[$i], "  ", $ys[$i], "  ", $cluster_obj->{label_cluster_best}->{$the_label}, "\n";
# }
}
