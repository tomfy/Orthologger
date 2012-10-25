#!/usr/bin/perl -w
use strict;

# This works as of Sept. 20, 2012. 
#use lib '/home/tomfy/Orthologger/lib';
use List::Util qw ( min max sum );
#use TlyUtil qw ( order_newick );
use CXGN::Phylo::Mrbayes;
# read in MrBayes .run1.t file
# remove branch lengths 
# put newicks in canonical form (so exactly one newick for each topology)
# store newicks in hash
# 
my $alignment_nex_filename = shift;		# e.g. fam98??.nex
my $burnin_frac = shift || 0.1;

my $mrb_obj = CXGN::Phylo::Mrbayes->new({'alignment_nex_filename' =>$alignment_nex_filename});
 print "# ", join("\n# ", @{$mrb_obj->get_id_species_string_array()}), "\n";
my $gen_param_hrefs = 
#load_params($alignment_nex_filename);
$mrb_obj->retrieve_param_samples(); # read in from *.run?.p file
my ($gen_ntopo_hrefs, $newick_number_map, $number_newick_map) = 
#  load_topologies($alignment_nex_filename);	#read in from * .run?.t file
$mrb_obj->retrieve_topology_samples();
# $topo_count is
my ($topo_count, $total_trees) = 
#  count_topologies($gen_ntopo_hrefs);
 $mrb_obj->count_topologies($gen_ntopo_hrefs);
my $distinct_newicks = scalar keys %$topo_count;
print "# Distinct topologies: $distinct_newicks\n";

my @sorted_tree_numbers = sort { sum(@{$topo_count->{$b}}) <=> sum(@{$topo_count->{$a}}) } keys %$topo_count; # tree numbers, sorted by occurrences (sum of runs).

my $number_id_map = 
#  load_number_id_map($alignment_nex_filename);
  $mrb_obj->retrieve_number_id_map();
print "# n trees (sorted): ", scalar @sorted_tree_numbers, "\n";
my  %number_rank_map = (); # number is n for nth distinct topology visited by chain
my %rank_number_map = (); # rank is 1 for topology with highest posterior prob. etc.
my $index = 1;
my $total_trees_so_far = 0;
my $sum_diff = 0;
my $L1_denom = 0; 
foreach my $tree_number (@sorted_tree_numbers) {
  my $treecount = sum @{$topo_count->{$tree_number}};
  my @treecounts = @{$topo_count->{$tree_number}};
  $L1_denom += $treecount;
  $sum_diff += max(@treecounts) - min(@treecounts);
  if ($treecount % 2 == 1) {
    $sum_diff--;
    $L1_denom--;
  }

  #  print "treecounts: ", join(" ", @treecounts), "\n";
  $total_trees_so_far += $treecount;
  my $newick_with_numbers = $number_newick_map->{$tree_number};
my $post_prob = $treecount/$total_trees;
  $number_rank_map{$tree_number} = $index;
  $rank_number_map{$index} = $tree_number;
  printf("%4i %4i %4i %4i %4i  %10.5g %10.5g  %4i  %s\n", 
	 $index, @treecounts, $treecount, $total_trees_so_far, 
	 $post_prob, $total_trees_so_far/$total_trees, $tree_number,
	 $newick_with_numbers);
  my $newick_with_ids = $mrb_obj->restore_ids_to_newick($newick_with_numbers, $number_id_map);
print $newick_with_ids, "\n";
  $index++;
}
#print "$sum_diff, $total_trees $L1_denom\n";
printf("# topology post. distrib. L1 difference between runs: %8.5f \n", 0.5*$sum_diff/$L1_denom);
print "# total tree hits: $total_trees \n";
print "# distinct newicks: $distinct_newicks \n";

my $j_run = 1;
foreach my $gen_ntopo (@$gen_ntopo_hrefs) {
  my $string = "# Run: $j_run\n# Gen Topo      LnL           TL           alpha         pinvar\n";
  my $trees_read_in = scalar keys %{$gen_ntopo};
print "# run: $j_run; trees read in: $trees_read_in.\n";
  my $n_burnin = int($burnin_frac * $trees_read_in);
  my @sorted_generations = sort {$a <=> $b} keys %{$gen_ntopo};
  # print "nburnin: $n_burnin  ", scalar @sorted_generations, "\n";
  foreach my $i_gen (@sorted_generations[$n_burnin..$trees_read_in-1]) {
    $string .= "$i_gen   " .
      $number_rank_map{$gen_ntopo->{$i_gen}} . "   " .
	$gen_param_hrefs->[$j_run-1]->{$i_gen} . "\n";
  }				# loop over gens
#  print "$string\n" if(1);
  open my $fhtp, ">run$j_run.tp";
  print $fhtp "$string\n";
  close $fhtp;
  $j_run++;
}				# loop over runs
# end of main
