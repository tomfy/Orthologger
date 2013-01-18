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
my $alignment_nex_filename = shift; # e.g. fam98??.nex
my $burnin_frac = shift || 0.1;
my $truncation_fraction = shift || 0.95;

my $mrb_obj = CXGN::Phylo::Mrbayes->new({'alignment_nex_filename' =>$alignment_nex_filename});
print "# ", join("\n# ", @{$mrb_obj->get_id_species_string_array()}), "\n";
my $gen_param_hrefs = 
  #load_params($alignment_nex_filename);
  $mrb_obj->retrieve_param_samples(); # read in from *.run?.p file
my ($generation_toponumber_hrefs, $newick_number_map, $number_newick_map) = 
  $mrb_obj->retrieve_topology_samples(); #read in from * .run?.t file
my $n_runs = scalar @$generation_toponumber_hrefs;
# $topo_count is
my ($topo_count, $total_trees) = 
  $mrb_obj->count_topologies($generation_toponumber_hrefs);
my $distinct_newicks = scalar keys %$topo_count;
print "# Distinct topologies: $distinct_newicks total trees: $total_trees\n";

# topo count: keys tree numbers, values array refs to arrays holding numbers of hits for each run.
my @sorted_tree_numbers = sort { sum(@{$topo_count->{$b}}) <=> sum(@{$topo_count->{$a}}) } keys %$topo_count; # tree numbers, sorted by occurrences (sum of runs).

my $number_id_map =
  #  load_number_id_map($alignment_nex_filename);
  $mrb_obj->retrieve_number_id_map();
print "# n trees (sorted): ", scalar @sorted_tree_numbers, "\n";
my  %number_rank_map = (); # number is n for nth distinct topology visited by chain
my %rank_number_map = (); # rank is 1 for topology with highest posterior prob. etc.
my $index = 1;
my $total_trees_so_far = 0;

if(1){
my $sum_diff = 0;
my $L1_denom = 0;
foreach my $tree_number (@sorted_tree_numbers) {
  my $treecount = sum @{$topo_count->{$tree_number}}; # number of hits for this topology, sum over runs.
  my @treecounts = @{$topo_count->{$tree_number}};
  # $L1_denom += $treecount;
  # $sum_diff += max(@treecounts) - min(@treecounts);
  # # if ($treecount % 2 == 1) { # ?? not sure this makes sense.
  # #   $sum_diff--;
  # #   $L1_denom--;
  # # }

  #  print "treecounts: ", join(" ", @treecounts), "\n";
  $total_trees_so_far += $treecount;
  my $newick_with_numbers = $number_newick_map->{$tree_number};
  my $post_prob = $treecount/$total_trees;
  $number_rank_map{$tree_number} = $index;
  $rank_number_map{$index} = $tree_number;
  printf("%4i ", $index);
  for (@treecounts) {
    printf("%4i ", $_);
  }
  printf("%4i %4i  %10.5g %10.5g  %4i  %s\n", 
	 $treecount, $total_trees_so_far, 
	 $post_prob, $total_trees_so_far/$total_trees, $tree_number
	 ,"" );
  # ,$newick_with_numbers);
  my $newick_with_ids = $mrb_obj->restore_ids_to_newick($newick_with_numbers, $number_id_map);
  # print $newick_with_ids, "\n";
  $index++;
}
}

print "# L1 (top ", 1-$truncation_fraction, " lumped):  ", L1_difference($topo_count, $truncation_fraction), "\n";;
print "# total tree hits: $total_trees \n";
print "# distinct newicks: $distinct_newicks \n";

my $j_run = 1;

foreach my $generation_toponumber (@$generation_toponumber_hrefs) {
  my $string = "# Run: $j_run\n# Gen Topo      LnL           TL           alpha         pinvar\n";
  my $trees_read_in = scalar keys %{$generation_toponumber};

  print "# run: $j_run; trees read in: $trees_read_in.\n";
  my $n_burnin = int($burnin_frac * $trees_read_in);
  my @sorted_generations = sort {$a <=> $b} keys %{$generation_toponumber};
  #my $n_runs = scalar @{$topo_count->{$generation_toponumber->{$sorted_generations[0]}}};
 
  # print "nburnin: $n_burnin  ", scalar @sorted_generations, "\n";
  foreach my $i_gen (@sorted_generations[$n_burnin..$trees_read_in-1]) {
    $string .= "$i_gen   " .
      $number_rank_map{$generation_toponumber->{$i_gen}} . "   " .
	$gen_param_hrefs->[$j_run-1]->{$i_gen} . "\n";
  }				# loop over gens
  #  print "$string\n" if(1);
  open my $fhtp, ">run$j_run.tp";
  print $fhtp "$string\n";
  close $fhtp;
  $j_run++;
}				# loop over runs
# end of main


sub L1_difference{
  my $toponumber_hitlist = shift;
  my $truncation_fraction = shift || 0.95;

  my @sorted_toponumbers = sort { sum(@{$toponumber_hitlist->{$b}}) <=> sum(@{$toponumber_hitlist->{$a}}) }
    keys %$toponumber_hitlist;	#
  my $total_hits = sum ( map(@{$toponumber_hitlist->{$_}}, @sorted_toponumbers) );
  my $run0_hits =  sum( map($toponumber_hitlist->{$_}->[0], @sorted_toponumbers));
  my $n_runs = scalar @{$toponumber_hitlist->{$sorted_toponumbers[0]}};
 
  my $L1_numerator = 0;
  my $L1_denominator = 0;
  my $cume_treecount = 0;
  my @cume_treecounts = ((0) x $n_runs);
  foreach my $tree_number (@sorted_toponumbers) { # loop over topologies
    my @treecounts = @{$toponumber_hitlist->{$tree_number}};
    @cume_treecounts = map { $cume_treecounts[$_] + $treecounts[$_] } 0..$#treecounts;
    my $treecount = sum(@treecounts); # number of hits for this topology, summed over runs.
    $cume_treecount += $treecount;

    @treecounts = sort { $b <=> $a } @treecounts;
    my $coeff = $n_runs - 1;
    for my $run_topo_hits (@treecounts) { # loop over runs
      my $sum_abs_diffs = $coeff*$run_topo_hits;
      $L1_numerator += $sum_abs_diffs;
      $L1_denominator += $run_topo_hits;
      $coeff -= 2;
    }				#print STDERR "\n";

    if ($cume_treecount >= $truncation_fraction*$total_hits) {
      my @other_treecounts = ();
      for (@cume_treecounts) {
	push @other_treecounts, $run0_hits - $_;
      }
      @other_treecounts = sort { $b <=> $a } @other_treecounts;
      $coeff =  $n_runs - 1;
      for my $run_other_hits (@other_treecounts) { # loop over runs
	my $sum_abs_diffs = $coeff*$run_other_hits;
	$L1_numerator += $sum_abs_diffs;
	$L1_denominator += $run_other_hits;
	$coeff -= 2;
      }
      last;
    }
  }
  print STDERR "# n runs: $n_runs \n";

  $L1_denominator *= ($n_runs > 1)? ($n_runs-1): 1;
  return $L1_numerator/$L1_denominator;
}
