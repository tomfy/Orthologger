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
my $alignment_nex_filename = shift;           # e.g. fam98??.nex
my $start_gen = shift || 0;
my $burnin_frac            = shift || 0.0;
my $lumped_tail_fraction   = shift || 0.05;
my $target_bin_weight         = shift || 0.02;


# my $zzz = equal_weight_bins( $target_bin_weight, @xxx );
# my $ebL1 = avg_L1_distance($zzz);
# print "ebL1: $ebL1 \n";


my $mrb_obj = CXGN::Phylo::Mrbayes->new(
    { 'alignment_nex_filename' => $alignment_nex_filename, 'burnin_fraction' => $burnin_frac } );
print "# ", join( "\n# ", @{ $mrb_obj->get_id_species_string_array() } ), "\n";
my $gen_param_hrefs =
  $mrb_obj->retrieve_param_samples();    # read in from *.run?.p file
my ( $generation_toponumber_hrefs, $newick_number_map, $number_newick_map ) =
  $mrb_obj->retrieve_topology_samples(undef, $start_gen);    #read in from * .run?.t file
my $n_runs = scalar @$generation_toponumber_hrefs;

# $topo_count is
my ( $topo_count, $total_trees ) =
  $mrb_obj->count_post_burn_in_topologies($generation_toponumber_hrefs);
my $distinct_newicks = scalar keys %$topo_count;
print "# Distinct topologies: $distinct_newicks total trees: $total_trees\n";

# topo count: keys tree numbers, values array refs to arrays holding numbers of hits for each run.
my @sorted_tree_numbers =
  sort { sum( @{ $topo_count->{$b} } ) <=> sum( @{ $topo_count->{$a} } ) }
  keys %$topo_count;    # tree numbers, sorted by occurrences (sum of runs).

my $number_id_map =
  #  load_number_id_map($alignment_nex_filename);
  $mrb_obj->retrieve_number_id_map();
print "# n trees (sorted): ", scalar @sorted_tree_numbers, "\n";
my %number_rank_map =
  ();                   # number is n for nth distinct topology visited by chain
my %rank_number_map =
  ();    # rank is 1 for topology with highest posterior prob. etc.
my $total_trees_so_far = 0;

my $index = 1;
foreach my $tree_number (@sorted_tree_numbers) {
    my $treecount =
      sum @{ $topo_count
          ->{$tree_number} }; # number of hits for this topology, sum over runs.
    my @treecounts = @{ $topo_count->{$tree_number} };
    $total_trees_so_far += $treecount;
    my $newick_with_numbers = $number_newick_map->{$tree_number};
    my $post_prob           = $treecount / $total_trees;
    $number_rank_map{$tree_number} = $index;
    $rank_number_map{$index}       = $tree_number;
    printf( "%4i ", $index );

    for (@treecounts) {
        printf( "%4i ", $_ );
    }
    printf(
        "%4i %4i  %10.5g %10.5g  %4i  %s\n",
        $treecount, $total_trees_so_far, $post_prob,
        $total_trees_so_far / $total_trees,
        $tree_number, ""
    );

    # ,$newick_with_numbers);
    my $newick_with_ids =
      $mrb_obj->restore_ids_to_newick( $newick_with_numbers, $number_id_map );

    # print $newick_with_ids, "\n";
    $index++;
}
# print "# L1alt (top $lumped_tail_fraction lumped):  ",
#   L1alt( $topo_count, $lumped_tail_fraction ), "\n";
my $avgL1_lrt = avg_L1_distance(CXGN::Phylo::Mrbayes::lump_right_tail( $topo_count, $lumped_tail_fraction ));
print "# L1 (top $lumped_tail_fraction lumped):  ", $avgL1_lrt, "\n";

my $y = rebin( $topo_count, $target_bin_weight );
print "# L1 (rebinned; fraction in each bin >= $target_bin_weight):  ",
  avg_L1_distance($y), "\n";

print "# total tree hits: $total_trees \n";
print "# distinct newicks: $distinct_newicks \n";

my $j_run = 1;
foreach my $generation_toponumber (@$generation_toponumber_hrefs) {
    my $string =
"# Run: $j_run\n# Gen Topo      LnL           TL           alpha         pinvar\n";
    my $trees_read_in = scalar keys %{$generation_toponumber};

    print "# run: $j_run; trees read in: $trees_read_in.\n";
    my $n_burnin = int( $burnin_frac * $trees_read_in );
    my @sorted_generations = sort { $a <=> $b } keys %{$generation_toponumber};

#my $n_runs = scalar @{$topo_count->{$generation_toponumber->{$sorted_generations[0]}}};

    # print "nburnin: $n_burnin  ", scalar @sorted_generations, "\n";
    foreach my $i_gen ( @sorted_generations[ $n_burnin .. $trees_read_in - 1 ] )
    {
        $string .=
            "$i_gen   "
          . $number_rank_map{ $generation_toponumber->{$i_gen} } . "   "
          . $gen_param_hrefs->[ $j_run - 1 ]->{$i_gen} . "\n";
    }    # loop over gens
         #  print "$string\n" if(1);
    open my $fhtp, ">run$j_run.tp";
    print $fhtp "$string\n";
    close $fhtp;
    $j_run++;
}    # loop over runs

# end of main

# doesn't do tail lumping or any sort of rebinning
sub avg_L1_distance {    # this one doesn't put the tail into an overflow bin
    my $label_weightslist = shift;

    my @labels     = keys %$label_weightslist;                               #
    my $total_hits = sum( map( @{ $label_weightslist->{$_} }, @labels ) );
    my $run0_hits  = sum( map( $label_weightslist->{$_}->[0], @labels ) );
    my $n_runs     = scalar @{ $label_weightslist->{ $labels[0] } };

    my $count_bins = 0;
    my $L1_numerator   = 0;
    my $L1_denominator = 0;
    my $cume_weight    = 0;
    my @cume_weights   = ( (0) x $n_runs );
    foreach my $label (@labels) {    # loop over topologies
        my @weights = @{ $label_weightslist->{$label} };
        @cume_weights = map { $cume_weights[$_] + $weights[$_] } 0 .. $#weights;
        my $weight =
          sum(@weights);   # number of hits for this topology, summed over runs.
        $cume_weight += $weight;

        @weights = sort { $b <=> $a } @weights;
        my $coeff = $n_runs - 1;
        for my $run_weight (@weights) {    # loop over runs
            my $sum_abs_diffs = $coeff * $run_weight;
            $L1_numerator   += $sum_abs_diffs;
            $L1_denominator += $run_weight;
            $coeff -= 2;
        }    #print STDERR "\n";
	$count_bins++;
    }
    print STDERR "# n runs: $n_runs  n bins: $count_bins  \n";

    $L1_denominator *= ( $n_runs > 1 ) ? ( $n_runs - 1 ) : 1;
    return $L1_numerator / $L1_denominator;
}

# sub right_tail_to_overflow {
#     my $label_weightslist = shift;
#     my $r_tail_weight     = shift || 0.05;
#     my %label_sumweights  = ();
#     while ( my ( $l, $ws ) = each %$label_weightslist ) {
#         $label_sumweights{$l} = sum(@$ws);
#     }
#     my @sorted_labels = sort

#    #   { sum(@{$label_weightslist->{$b}}) <=> sum(@{$label_weightslist->{$a}}) }
#     {
#         $label_sumweights{$b} <=> $label_sumweights{$a}
#       }
#       keys %$label_weightslist;    #
#     my $total_hits =
#       sum( map( @{ $label_weightslist->{$_} }, @sorted_labels ) );
#     my $run0_hits = sum( map( $label_weightslist->{$_}->[0], @sorted_labels ) );
#     my $n_runs = scalar @{ $label_weightslist->{ $sorted_labels[0] } };

#     my $result       = {};
#     my $cume_weight  = 0;
#     my @cume_weights = ( (0) x $n_runs );
#     foreach my $label (@sorted_labels) {    # loop over categories
#         my @weights = @{ $label_weightslist->{$label} };
#         @cume_weights = map { $cume_weights[$_] + $weights[$_] } 0 .. $#weights;
#         my $weight =
#           sum(@weights); # number of hits for this categories, summed over runs.
#         $cume_weight += $weight;
#         $result->{$label} = \@weights;
#         if ( $cume_weight >= ( 1 - $r_tail_weight ) * $total_hits ) {
#             my @other_weights = ();
#             for (@cume_weights) {
#                 push @other_weights, $run0_hits - $_;
#             }
#             $result->{ $label + 1 } = \@other_weights;
#             last;
#         }
#     }
#     print STDERR "# n runs: $n_runs \n";
#     return $result;
# }

sub rebin # input here is already binned
{ # the idea here is to make each bin have approx. same fraction of total weight (1% is default)
    my $label_weightslist = shift;
    my $target_bin_weight = shift || 0.01;

    my %label_sumweights = ();
    while ( my ( $l, $ws ) = each %$label_weightslist ) {
        $label_sumweights{$l} = sum(@$ws);
    }
    my @sorted_labels = sort {
        $label_sumweights{$a} <=> $label_sumweights{$b}
      }                            # sort by weight; small to large
      keys %$label_weightslist;    #

    my $total_hits =
      sum( map( @{ $label_weightslist->{$_} }, @sorted_labels ) );
    my $run0_hits = sum( map( $label_weightslist->{$_}->[0], @sorted_labels ) );
    my $n_runs = scalar @{ $label_weightslist->{ $sorted_labels[0] } };

    my $result       = {};
    my $cume_weight  = 0;
    my @cume_weights = ( (0) x $n_runs );
    foreach my $label (@sorted_labels) {    # loop over categories
        my @weights = @{ $label_weightslist->{$label} };
        @cume_weights = map { $cume_weights[$_] + $weights[$_] } 0 .. $#weights;
        my $weight =
          sum(@weights); # number of hits for this categories, summed over runs.
        $cume_weight += $weight;

        #  $result->{$label} = \@weights;
        # print "$weight, $cume_weight, $total_hits, $target_bin_weight \n";
        if ( $cume_weight >= $target_bin_weight * $total_hits ) {
            my @copy = @cume_weights;
            $result->{$label} = \@copy;
 #           print join( " ", @cume_weights ), "  ", sum(@copy), "\n";
            $cume_weight = 0;
            @cume_weights = ( (0) x $n_runs );
        }
    }
    print "number of bins ", scalar keys %$result, "\n";

#  print STDERR "# n runs: $n_runs \n\n\n";
# for(keys %$result){
#  print "$_ ", join("; ", @{$result->{$_}}), "  ", sum(@{$result->{$_}}), "\n";
# }
    return $result;
}



sub equal_weight_bins {    #input here is not yet binned, just a few sets of data, each
# just an array of N data points (numbers).
    my $min_bin_fraction = shift || 0.01;
    my @data_sets =
      @_;    # each element is array ref storing data points (numbers).

   my $result = {};
    my @data_set_maxes = ();
    for (@data_sets) {
        push @data_set_maxes, max(@$_);
        my @sorted_data_set = sort { $a <=> $b } @$_;
      #  print STDERR "A: ", join( ", ", @sorted_data_set ), "\n";

        #    push @sorted_data_set, undef;
        $_ = \@sorted_data_set;

        #print STDERR "B: ", join(", ", @sorted_data_set), "\n";
    }
    my $max_data = max(@data_set_maxes);
    my $big      = $max_data + 1e100;
    for (@data_sets) {
        push @$_, $big;
    }
 #   print STDERR join( ";", @data_set_maxes ), "  $max_data \n";
    my $n_data_sets     = scalar @data_sets;               #
    my @n_points        = map( (scalar @$_-1), @data_sets );
    my $n_points_in_set = min(@n_points);
    warn "Different numbers of data points in different runs: ",
      join( ", ", @n_points ), "\n"
      if ( min(@n_points) != max(@n_points) );
    my $n_total_points = $n_data_sets * $n_points_in_set;

    my $i                    = 1;
    my $points_binned_so_far = 0;
    my $bin_number           = 0;
    my @next_point_indices   = ( (0) x $n_data_sets );
    my @counts_this_bin      = ( (0) x $n_data_sets );
    my $total_this_bin       = 0;
    while (1) {
        my $desired_cumulative_points =
          int( $i * $min_bin_fraction * $n_total_points + 0.5 );
        my @next_points = ();
        while ( my ( $i, $points ) = each @data_sets ) {
            $next_points[$i] = $points->[ $next_point_indices[$i] ];
        }
        my ( $i_min, $min ) = ( 0, $next_points[0] );
        while ( my ( $i, $v ) = each(@next_points) ) {
     #       print STDERR "i: v:  $i  ", ( defined $v ) ? $v : 'undef', " \n";
            if ( defined $v and $v <= $min ) {
                $min   = $v;
                $i_min = $i;
           #     print STDERR "imin min: $i_min  $min \n";
            }
        }
    #    print STDERR "\n";
        last if ( $min > $max_data );
        $counts_this_bin[$i_min]++;
        $points_binned_so_far++;
        $total_this_bin++;
        $next_point_indices[$i_min]++;
#	print STDERR "YYYYYYY: $points_binned_so_far,  $desired_cumulative_points\n";
        if ( $points_binned_so_far >= $desired_cumulative_points ) {
            my @copy = @counts_this_bin;
            $result->{$bin_number} = \@copy;
            @counts_this_bin = ( (0) x $n_data_sets );
            $total_this_bin = 0;
            $bin_number++;
	    $i++;
        }
      }
    # for ( keys %$result ) {
    #     print "$_     ", join( ", ", @{ $result->{$_} } ), "\n";
    # }
    return $result;

}
##########################################################
# not tested
sub right_tail_to_overflow_alt
{    # arg is list of hashrefs with label/weight pairs
    my @histograms =
      @{ (shift)
      }; # elements are hashrefs; keys are category labels, values are weights (probabilities, counts,...)
    my $r_tail_weight = shift || 0.05;
    my @histogram_weight_sums = map ( sum( values %$_ ), @histograms );

    # for(@histograms){
    #   push @sums, sum(values %$_);
    # }

    for (@histogram_weight_sums) {
        if ( $_ != $histogram_weight_sums[0] ) {
            warn "Histogram sums not all same: ",
              join( "; ", @histogram_weight_sums ), "\n";
            last;
        }
    }
    my $total_weight = sum(@histogram_weight_sums);
    my $n_histograms = scalar @histograms;

    my $sum_of_histograms;
    for my $a_histogram (@histograms) {
        while ( my ( $label, $weight ) = each %$a_histogram ) {
            $sum_of_histograms->{$label} += $weight;
        }
    }

    my @sorted_labels =
      sort { $a <=> $b } keys %$sum_of_histograms;    # smallest label at [0]
    my $cumulative_weight = 0;
    my $in_tail           = 0;
    my $rightmost_regular_bin_label;
    my @overflows = ( (0) x $n_histograms );
    while ( my ( $i, $label ) = each @sorted_labels ) {
        $cumulative_weight += $sum_of_histograms->{$label};
        if ($in_tail) {
            while ( my ( $i, $value ) = each @overflows )
            {    # move bin contents into overflow bin
                $overflows[$i] += $histograms[$i]->{$label};
                delete $histograms[$i]->{$label};
            }
        }
        else {
            $rightmost_regular_bin_label = $label;
            if ( $cumulative_weight > $r_tail_weight * $total_weight ) {
                $in_tail = 1
                  ; # the rest of the bins will eliminated, their contents put into overflow.
            }
        }
    }
    my $overflow_label = $rightmost_regular_bin_label + 1;
    while ( my ( $i, $histogram ) = each @histograms ) {
        $histogram->{$overflow_label} = $overflows[$i];
    }
    for my $histogram (@histograms) {
        for ( sort { $a <=> $b } keys %$histogram ) {
            print STDERR "$_  ", $histogram->{$_}, "\n";
        }
    }
    return \@histograms;
}

#  obsolete -

# tail lumping done in subroutine
sub L1alt {
    my $label_weightslist = shift;
    my $r_tail_weight = shift || 0.05;

    my @sorted_labels = sort {
        sum( @{ $label_weightslist->{$b} } )
          <=> sum( @{ $label_weightslist->{$a} } )
      }
      keys %$label_weightslist;    #
    my $total_hits =
      sum( map( @{ $label_weightslist->{$_} }, @sorted_labels ) );
    my $run0_hits = sum( map( $label_weightslist->{$_}->[0], @sorted_labels ) );
    my $n_runs = scalar @{ $label_weightslist->{ $sorted_labels[0] } };

    my $L1_numerator   = 0;
    my $L1_denominator = 0;
    my $cume_weight    = 0;
    my @cume_weights   = ( (0) x $n_runs );
my $count_bins = 0;
    foreach my $label (@sorted_labels) {    # loop over topologies
        my @weights = @{ $label_weightslist->{$label} };
        @cume_weights = map { $cume_weights[$_] + $weights[$_] } 0 .. $#weights;
        my $weight =
          sum(@weights);   # number of hits for this topology, summed over runs.
        $cume_weight += $weight;

        @weights = sort { $b <=> $a } @weights;
        my $coeff = $n_runs - 1;
        for my $run_weight (@weights) {    # loop over runs
            my $sum_abs_diffs = $coeff * $run_weight;
            $L1_numerator   += $sum_abs_diffs;
            $L1_denominator += $run_weight;
            $coeff -= 2;
	  
        }    #print STDERR "\n";
	$count_bins++;
        if ( $cume_weight >= ( 1 - $r_tail_weight ) * $total_hits ) {
            my @other_weights = ();
            for (@cume_weights) {
                push @other_weights, $run0_hits - $_;
            }
            @other_weights = sort { $b <=> $a } @other_weights;
            $coeff = $n_runs - 1;
            for my $run_other_hits (@other_weights) {    # loop over runs
                my $sum_abs_diffs = $coeff * $run_other_hits;
                $L1_numerator   += $sum_abs_diffs;
                $L1_denominator += $run_other_hits;
                $coeff -= 2;
            }
	    $count_bins++;
            last;
        }
    }
    print STDERR "# n runs, n bins: $n_runs $count_bins \n";

    $L1_denominator *= ( $n_runs > 1 ) ? ( $n_runs - 1 ) : 1;
    return $L1_numerator / $L1_denominator;
}
