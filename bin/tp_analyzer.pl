#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );
# read in MrBayes .run1.t file
# remove branch lengths 
# store newicks in hash

my $pattern = shift;		# e.g. fam98??.nex
my $pattern1 = "$pattern.run?";

my $burnin_frac = shift || 0.1;

my $gen_param_hrefs = load_params($pattern);

my ($gen_ntopo_hrefs, $newick_number_map, $number_newick_map) = 
    load_topologies($pattern);

my ($topo_count, $total_trees) = count_topologies($gen_ntopo_hrefs);

my $distinct_newicks = scalar keys %$topo_count;
print "Distinct topologies: $distinct_newicks\n";

my @trees = sort { sum(@{$topo_count->{$b}}) <=> sum(@{$topo_count->{$a}}) } keys %$topo_count;

print "n trees (sorted): ", scalar @trees, "\n";
my  %number_rank_map = ();
my %rank_number_map = ();
my $index = 1;
my $total_trees_so_far = 0;
foreach (@trees) {
  my $treecount = sum @{$topo_count->{$_}};
  my @treecounts = @{$topo_count->{$_}};
#  print "treecounts: ", join(" ", @treecounts), "\n";
  $total_trees_so_far += $treecount;
  my $post_prob = $treecount/$total_trees;
$number_rank_map{$_} = $index;
$rank_number_map{$index} = $_;
  printf("%4i %4i %4i %4i %4i  %10.5g %10.5g  %4i  %s\n", 
	 $index, @treecounts, $treecount, $total_trees_so_far, 
	 $post_prob, $total_trees_so_far/$total_trees, $_,
	 $number_newick_map->{$_});
  $index++;
}

print "total tree hits: $total_trees \n";
print "distinct newicks: $distinct_newicks \n";
#

my $j_run = 1;
foreach my $gen_ntopo (@$gen_ntopo_hrefs) {
  my $string = "Run: $j_run\n";
  my $trees_read_in = scalar keys %{$gen_ntopo};
print "run, trees read in $j_run  $trees_read_in \n";
  my $n_burnin = int($burnin_frac * $trees_read_in);
  my @sorted_generations = sort {$a <=> $b} keys %{$gen_ntopo};
  # print "nburnin: $n_burnin  ", scalar @sorted_generations, "\n";
  foreach my $i_gen (@sorted_generations[$n_burnin..$trees_read_in-1]) {
    $string .= "$i_gen " .
      $number_rank_map{$gen_ntopo->{$i_gen}} . " " .
	$gen_param_hrefs->[$j_run-1]->{$i_gen} . "\n";
  }				# loop over gens
#  print "$string\n" if(1);
  open my $fhtp, ">run$j_run.tp";
  print $fhtp "$string\n";
  close $fhtp;
  $j_run++;
}				# loop over runs
# end of main

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


sub load_params{
  # read data from  run?.p file
  # store in 
  my $pattern = shift;		# e.g. fam9877.nex
  my $p_files = `ls $pattern.run?.p`;
  my @p_infiles = split(" ", $p_files);
  my $n_runs_p = scalar @p_infiles;

  # the following has one elem for each run, and it is
  # a hash ref, with generation numbers as keys, 
  # parameter strings (logL, TL, alpha ...) as values
  my @gen_param_hashes = ();
  foreach (1..$n_runs_p) {
    push @gen_param_hashes, {};
  }
  #my $p_run = 1;
  foreach my $i_run_p (1..$n_runs_p) {
    my $p_file = "$pattern.run$i_run_p.p";
    open  my $fhp, "<$p_file";

    while (my $line = <$fhp>) {
      chomp $line;
      next unless($line =~ /^\s*(\d+)/);
      #  print "$line \n";
      my @cols = split(" ", $line);
      my $generations = shift @cols;
      my $param_string = join("  ", @cols);
      $gen_param_hashes[$i_run_p-1]->{$generations} = $param_string;
    }
    $i_run_p++;
  }
  return \@gen_param_hashes;
}
# end of reading in parameter values

sub load_topologies{
my $pattern = shift; # e.g. fam987?.nex

 my $t_files = `ls $pattern.run?.t`;
  my @t_infiles = split(" ", $t_files);
  my $n_runs = scalar @t_infiles;
  my @gen_ntopo_hashes = ();
  foreach (1..$n_runs) {
    push @gen_ntopo_hashes, {};
  }
my %newick_number_map = ();
my %number_newick_map = ();
  my $topology_count = 0;
  foreach my $i_run (1..$n_runs) {
    my $t_infile = "$pattern.run$i_run.t";
    open my $fh, "<$t_infile";

    # read trees in, remove branch lengths, store in array

    while (my $line = <$fh>) {
      chomp $line;
      if ($line =~ s/tree gen[.](\d+) = .{4}\s+//) {
	my $newick = $line;
	my $generation = $1;
	$newick =~ s/:[0-9e\-.]*(,|[)])/$1/g; # remove branch lengths
	#print "[$newick]\n";
	$newick =~ s/^\s+//; 
	$newick =~ s/;\s*//;
	$newick = order_newick($newick);
	#	print $newick, "\n";
	#		exit;
	if (!exists $newick_number_map{$newick}) {
	  $topology_count++;
	  $newick_number_map{$newick} = $topology_count; # 1,2,3,...
	  $number_newick_map{$topology_count} = $newick;
	}
	$gen_ntopo_hashes[$i_run-1]->{$generation} = $newick_number_map{$newick};
      }
    } # now $gen_ntopo_hashes[$i_run] is hash ref with generations as keys, and topology numbers as values.
  }
return (\@gen_ntopo_hashes, \%newick_number_map, \%number_newick_map);
}

sub count_topologies{
my $gen_ntopo_hrefs = shift;
my %topo_count = ();
my $total_trees = 0;
foreach my $i_run (1..scalar @$gen_ntopo_hrefs) {
  my $gen_ntopo = $gen_ntopo_hrefs->[$i_run-1];
  my $trees_read_in = scalar keys %{$gen_ntopo};
  print "Run: $i_run. Trees read in: $trees_read_in\n";
  # store trees from array in hash, skipping burn-in
  my $n_burnin = int($burnin_frac * $trees_read_in);
#  print "trees read in: $trees_read_in. Post burn-in: ", $trees_read_in - $n_burnin, "\n";
  my @sorted_generations = sort {$a <=> $b} keys %{$gen_ntopo};
  foreach my $i_gen (@sorted_generations[$n_burnin..$trees_read_in-1]) {
    my $topo_number = $gen_ntopo->{$i_gen};
    if (!exists $topo_count{$topo_number}) {
      my @zeroes = ((0) x scalar @$gen_ntopo_hrefs);
      $topo_count{$topo_number} = \@zeroes;
    }
    $topo_count{$topo_number}->[$i_run-1]++;
    $total_trees++;
  }
  $i_run++;
}
return (\%topo_count, $total_trees);
}
