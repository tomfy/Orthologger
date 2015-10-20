#!/usr/bin/perl -w
use strict;
use List::Util qw( min max sum );
use Graph;

# read in blast output (abc format) and group sequences into families
# basic idea: The set of matches to each query in the abc file defines a 'family'
# however there may be a lot of overlap between these families; so we want to
# join together highly overlapping families to make fewer families.
# two query sequences are considered connected if a fraction >= $min_overlap_fraction 
# of the matches in the smaller of the two sets of matches are also present in the
# larger match set.
# Then construct the graph with edges joining connected pairs, and find the connected subgraphs;
# Each connected subgraph has a set of nodes which are query sequences; these plus the union of
# all of their sets of matches make up the family.

my $min_overlap_fraction = shift || 0.9;
print "# min overlap fraction: $min_overlap_fraction \n";

my %query_match = (); # key: query id, value: hashref of match_id:e-value pairs
my %match_query = (); # key: match id, value: hashref of query_id:e-value pairs 

# Read in abc file and store in hashes.
my $lines_read = 0;
while (my $abc_line = <>) {
  $lines_read++;
  my ($id1, $id2, $e_value) = split(" ", $abc_line);
  if (! exists $query_match{$id1}) {
    $query_match{$id1} = {$id2 => $e_value};
  } else {
    $query_match{$id1}->{$id2} = $e_value;
  }
  if (! exists $match_query{$id2}) {
    $match_query{$id2} = {$id1 => $e_value}; # $id2 is new key, set value to hashref with 1 key:value pair
  } else {
    $match_query{$id2}->{$id1} = $e_value; # 
  }
  if ($lines_read == 1000) {
     print STDERR "Done storing  ", scalar keys %query_match, " queries.\n"; # if(($lines_read % 10000) == 0);
    $lines_read = 0;
  }
}
printf STDERR "Done storing abc info in hashes. \n";
# Create graph reflecting connectedness of queries
my $count_vertices = 0;
my $G = Graph::Undirected->new();
my $edge_count = 0;
my %connected_pairs = (); # keys "$id_a $id_b" where the ids are in lex. order in each pair, value = 1
for my $id1 (keys %query_match) {
  $G->add_vertex($id1);
  $count_vertices++;
  # only consider for id2 other queries which are matches of id1.
  my $conn_to_id1 = connected_queries($id1, \%query_match, \%match_query, $min_overlap_fraction, \%connected_pairs);
print STDERR "adding edges for queries connected to query $id1 \n";
  for my $id2 (@$conn_to_id1) {
    $edge_count++;
    $G->add_edge($id1, $id2);
  }
}
print "# conn pairs: ", scalar keys %connected_pairs, " sum of values:  ", sum(values %connected_pairs), "\n";

# Get the connected components; each of which is a set of queries which go in the same family.
my $sum_sizes = 0;
print STDERR "Now get connected components of graph.\n";
my @connected_components = $G->connected_components();
print STDERR "Done getting connected components. \n";
print "# Number of connected components: ", scalar @connected_components, "\n";
for my $ccomponent (@connected_components) { # loop over connected components (i.e. families)
  my %id_in_component__q_match_count = (); # keys are ids in family (not just core ids), values are count (of queries that the id matches)
  my %match_minev = (); # keys are ids in family (not just core), values are min e-values over core ids.
  my $biggest_q_match_set_size = -1;
  my $ccsize = scalar @$ccomponent; # @$ccomponent is array of query ids
  $sum_sizes += $ccsize; # total number of queries which have been dealt with so far.
  my $sum_id1_matches = 0;
  for my $core_id (@$ccomponent) { # for each core id ...
    my $q_match_set_size = scalar keys %{$query_match{$core_id}};
    $sum_id1_matches += $q_match_set_size;
    $biggest_q_match_set_size = max($q_match_set_size, $biggest_q_match_set_size);
    for my $match_id (keys %{$query_match{$core_id}}) { # increment the count of each id which is a match to it:\
      my $match_evalue = $query_match{$core_id}->{$match_id};
      $id_in_component__q_match_count{$match_id}++;
      $match_minev{$match_id} = (exists $match_minev{$match_id})? 
	min($match_minev{$match_id}, $match_evalue) :
	  $match_evalue;
    }
  }
  my @vals = values %id_in_component__q_match_count;
  my $mult_count = histogram(\@vals);
  my @sskeys = sort {$b <=> $a} keys %$mult_count;


  my $fam_size = scalar @vals;
  my ($overlap90_count, $overlap80_count, $overlap50_count) = (0, 0);
  my $count_all = 0;
  for (@sskeys) {
  #  print "$_: ", $mult_count->{$_}, "; ";
    if ($_ >= 0.9*$ccsize) { # number of seqs which match at least 90% of the core ids.
      $overlap90_count += $mult_count->{$_};
    }
  if ($_ >= 0.8*$ccsize) { # number of seqs which match at least 80% of the core ids.
      $overlap80_count += $mult_count->{$_};
    }
    if ($_ >= 0.5*$ccsize) { # number of seqs which match at least 50% of the core ids.
      $overlap50_count += $mult_count->{$_};
    }
    $count_all += $mult_count->{$_};
  }# print "\n";
  my $union_size = scalar keys %match_minev;
  
  my $cc_size_all =  scalar keys %id_in_component__q_match_count;
  print "[ ", join(" ", sort @$ccomponent), " ]   "; # the queries of the nodes in connected component 
  printf("%4i %5i %4.2f %4.2f %4.2f %4.2f\n", $ccsize, $cc_size_all, $sum_id1_matches/($ccsize*$cc_size_all),
	 $overlap90_count/$count_all, 
#	 $overlap80_count/$count_all, 
	 $overlap50_count/$count_all,
	  $biggest_q_match_set_size/$union_size); # if($ccsize > 1);
  for (keys %match_minev){
    print "    $_  ", $match_minev{$_}, "\n";
  }print "\n";
} # end of loop over connected components.
print "# Number of connected components: ", scalar @connected_components, "\n";
print "# vertex count: $count_vertices  edge count: $edge_count  sum of connected component sizes: $sum_sizes \n";

sub connected_queries{ # given a query id, return a array ref of other queries connected to it.
  # but only the ones which weren't previously found and stored in $connected_pairs 
  my $query = shift;
  my $q_m_href = shift;
  my $m_q_href = shift;
  my $min_overlap_fraction = shift;
  my $connected_pairs = shift;

  my @connected_to_q = ();
  for my $match ( keys %{$q_m_href->{$query}} ) {
    next if ($query eq $match); # don't need to put matches with self into graph.
    if (exists $q_m_href->{$match}) { # $match  is a key of %$q_m_href, i.e. is itself a query id.
      my $pair_id = ($query le $match)? "$query $match" : "$match $query";
      next if (exists $connected_pairs->{$pair_id}); # { # this pair already known to be connected.
      my $other_query = $match;
      if (query_pair_connected($query, $other_query, $q_m_href, $m_q_href, $min_overlap_fraction)) {
	push @connected_to_q, $other_query;
	$connected_pairs->{$pair_id}++;
      }
    }
  }				# end of loop
  return \@connected_to_q;
}

sub query_pair_connected{ # does not use e-value - just uses presence/absence in q_m_href and m_q_href.
  my $qid1 = shift;
  my $qid2 = shift;
  my $q_m_href = shift;
  my $m_q_href = shift;
  my $min_overlap_fraction = shift || 0.9;

  return 0 if(!exists $q_m_href->{$qid1}->{$qid2}  or  !exists $q_m_href->{$qid2}->{$qid1});

  my $size1 = scalar keys  %{$q_m_href->{$qid1}}; # size of set of matches to $qid1
  my $size2 = scalar keys  %{$q_m_href->{$qid2}};
  my $Ok = 0;
  if ($size1 <= $size2) { # overlap defined as fraction >= $min_overlap_fraction
# of smaller set is in intersection with larger set.
    $Ok = overlap($qid1, $qid2, $q_m_href, $m_q_href, $min_overlap_fraction);
  } else {
    $Ok = overlap($qid2, $qid1, $q_m_href, $m_q_href, $min_overlap_fraction);
  }
  return $Ok;
}

# Counts how many are matches to both qid1 and qid2, and how many are matches to qid1, but not qid2.
sub overlap{ # $qid1 is the id of the query with the smaller set of matches.
  my $qid1 = shift;
  my $qid2 = shift;
  my $q_m_href = shift;
  my $m_q_href = shift;
  my $min_overlap_fraction = shift;

  my $size1 = scalar keys %{$q_m_href->{$qid1}};
  my $max_1not2 = int((1 - $min_overlap_fraction)*$size1);
  my ($count_both, $count_1not2) = (0, 0);
  while (my ($mid1, $ev1) = each %{$q_m_href->{$qid1}}) {
    if (exists $m_q_href->{$mid1}->{$qid2}) { # i.e. $qid2 has $mid1 as a match
      $count_both++;		# both ev1, ev2 are < threshold
    } else {
      $count_1not2++;
      return 0 if($count_1not2 > $max_1not2);
    }
  }
  return 1;
}

sub histogram{			# histogram an array of integers.
  my $aref = shift;
  my %hist = ();
  for (@$aref) {
    $hist{$_}++;
  }
  return \%hist;
}

# sub big_and_union_sizes{
# my $connected_component = shift;
# my @cc = @$connected_component; # array of vertices
# for(@cc){



sub count_cc_edges{
  my $conncomp = shift;		# ref to array of vertices
  my $conn_pairs = shift;	# keys: "id1 id2"
  my $edge_count = 0;	 # counts the edges in the connected component
  my @cc = @$conncomp;
  my $ccsize = scalar @cc;
  my $max_edges = $ccsize*($ccsize - 1)/2;
  for (my $i = 0; $i < $ccsize; $i++) {
    my $id1 = $cc[$i];
    for (my $j = $i+1; $j < $ccsize; $j++) {
      my $id2 = $cc[$j];
      my $idpair = ($id1 le $id2)? "$id1 $id2" : "$id2 $id1";
      $edge_count++ if(exists $conn_pairs->{$idpair});
    }
  }
  return ($edge_count, $max_edges);
}

