#!/usr/bin/perl -w
use strict;
use List::Util qw( min max sum );
use Graph;

# read in file with amangout format
# query id in 1st column, and col --
# has [id1,id2,id3,...idn] a string with the ids in the AMangiosperm clade
# store these, then do pairwise comparisons, make graph with node for each
# query id, nodes connected if corresponding AMang clade id sets are
# highly overlapping.
my $biconnected = shift || 0;
my $verbose = shift || 0;
my $min_overlap_fraction = shift || 0.25;
print "# min overlap fraction: $min_overlap_fraction \n";

my %qid_clamangids = ();
my @qids = ();
while (<>) {
   next if(/^\s*#/);
   my @cols = split(" ", $_);
   my ($qid, $ids) = @cols[0, 6];
   push @qids, $qid;
   if ( $ids =~ /(\[.*\])\s*(\d+)\s+\d/) {
      $ids = $1;
   }
   $ids =~ s/(\[|\])//g;
   my @clade_ids = split(",", $ids);
   # my %clade_id_hash = ();
   # for (@clade_ids) {
   #    $clade_id_hash{$_} = 1;
   # }
   if (!exists $qid_clamangids{$qid}) {
      $qid_clamangids{$qid} = {};
   }
   for my $an_id (@clade_ids) {
      $qid_clamangids{$qid}->{$an_id} = 1;
   }
   # $qid_clamangids{$qid} = \%clade_id_hash;
   # print "$qid  ", join("; ", @clade_ids), "\n";
}
# done reading from file

my $G = Graph::Undirected->new();
# find the overlaps
my %q1__q2_overlap = ();
for my $qid1 (@qids) {
   $q1__q2_overlap{$qid1} = {};
   for my $qid2 (@qids) {
      my ($n_1only, $n_both, $n_2only) = venn2clade($qid1, $qid2, \%qid_clamangids);
      my $f_overlap = $n_both/($n_both + min($n_1only, $n_2only));
      if ($f_overlap > $min_overlap_fraction) {
         # in case there are multiple families with same qid (e.g. made using different alignment methods, mafft, muscle )
         # the overlap between two qids is the max overlap among their 
         $q1__q2_overlap{$qid1}->{$qid2} = (!exists $q1__q2_overlap{$qid1}->{$qid2})? $f_overlap : max($q1__q2_overlap{$qid1}->{$qid2}, $f_overlap);
         printf("%30s %30s  %4i  %4i %4i %4i  %4.3f \n", $qid1, $qid2, $f_overlap, $n_1only, $n_both, $n_2only, $q1__q2_overlap{$qid1}->{$qid2}) if($f_overlap > $min_overlap_fraction);
      }
   }
   my $qid2_overlap = $q1__q2_overlap{$qid1};
   my @sqid2s = sort {$qid2_overlap->{$b} <=> $qid2_overlap->{$a}} keys %$qid2_overlap;
   if ($verbose) {
      for (@sqid2s) {
         print "# $qid1  $_  ", $qid2_overlap->{$_}, "\n";
         last if(!$verbose);
      }                         # print "\n";
   }
   $G->add_vertex($qid1);
}
my $singleton_count = 0;
while (my ($q1, $q2s) = each %q1__q2_overlap) {
 #  print join(";;", keys %$q2s), "\n";
   $singleton_count++ if(scalar keys %$q2s == 1);
}

for my $qid1 (@qids) {
   for my $qid2 (keys %{$q1__q2_overlap{$qid1}}) {
      $G->add_edge($qid1, $qid2) if($qid1 ne $qid2);
      # $G->add_edge($qid1, $qid2);
   }
}
my @connected_components = $biconnected? $G->biconnected_components(): $G->connected_components();

# my %cc_clamangids = ();
my $sum_cc_sizes = 0;
my $cc_number = 1;
for my $ccomponent (@connected_components) { # loop over connected components (i.e. families)
   # $ccomponent is an array ref of vertices
   $sum_cc_sizes += scalar @$ccomponent;
   my $min_degree = 1000000;
   my $max_degree = -1;
   my @sverts = sort @$ccomponent;
   my $union_clamangids = {};
   for my $v (@sverts) {       #(@$ccomponent) {
      #    print "$v ", $G->degree($v), "\n",
      $min_degree = min($min_degree, $G->degree($v));
      $max_degree = max($max_degree, $G->degree($v));
      my $clamangids = $qid_clamangids{$v};
      for(keys %$clamangids){
         $union_clamangids->{$_}++;
      }
   }                            # print "\n";
 #  $cc_clamangids{$cc_number} = $union_clamangids;

   print "conn. component size: ", scalar @$ccomponent, "  min, max degree: $min_degree $max_degree  ", join(",", @sverts), " ",
     scalar keys %$union_clamangids, " [", join(",", keys %$union_clamangids), "]",
       "\n";
   # if($min_degree < scalar @$ccomponent-1  or  $max_degree >  scalar @$ccomponent-1){
   #    print join(" ", @$ccomponent), "\n";
   # }
   $cc_number++;
} # loop over connected components
my $n_conn_comps = scalar @connected_components;
if ($biconnected) {
   $n_conn_comps += $singleton_count;
   $sum_cc_sizes += $singleton_count;
}
print "# n conn. comps: $n_conn_comps   sum cc sizes: $sum_cc_sizes  singletons:  $singleton_count \n";

# the end

sub venn2clade{
   my $qid1 = shift;
   my $qid2 = shift;
   my $clade_id_href  = shift;
   my $clade1_ids = $clade_id_href->{$qid1};
   my $clade2_ids = $clade_id_href->{$qid2};
   my $clade1_size = scalar keys %$clade1_ids;
   my $clade2_size = scalar keys %$clade2_ids;
   my $n_intersect = 0; # counts number of ids in intersection of clade1, clade2
   if ($clade1_size <= $clade2_size) {
      for my $id1 (keys %$clade1_ids) {
         $n_intersect++ if(exists $clade2_ids->{$id1});
      }
   } else {
      for my $id2 (keys %$clade2_ids) {
         $n_intersect++ if(exists $clade1_ids->{$id2});
      }
   }
   my ($only1, $only2) = ($clade1_size - $n_intersect, $clade2_size - $n_intersect);
   return ($only1, $n_intersect, $only2);
}




# my %query_qAMids = (); # key: query id, value: hashref of id:1 pairs, one for each id in qAM subtree
 

# # Read in apuwids file and store in hashes.
# my $lines_read = 0;
# while (my $apuwids_line = <>) {
#    $lines_read++;
#    next if($apuwids_line =~ /^\s*$/); # skip whitespace-only lines
#    next if($apuwids_line =~ /^\s*#/); # skip comments
#    $apuwids_line =~ /^\s*(\S+).*\[(.*)\]/;
#    my $qid = $1;
#    $query_qAMids{$qid} = {};
#    my @ids = split(" ", $2); 
#    for (@ids) {
#       $query_qAMids{$qid}->{$_} = 1;
#    }
# }
# printf STDERR "Done storing apuwids info in hashes. \n";


# # Create graph reflecting connectedness of queries
# my $count_vertices = 0;
# my $G = Graph::Undirected->new();
# my $edge_count = 0;
# my %connected_pairs = (); # keys "$id_a $id_b" where the ids are in lex. order in each pair, value = 1
# #my $count = 0;
# for my $id1 (keys %query_qAMids) {
#    $G->add_vertex($id1);
#    $count_vertices++;
#    # only consider for id2 other queries which are matches of id1.
#    my $conn_to_id1 = connected_queries($id1, \%query_qAMids, $seqid_species, $min_overlap_fraction, \%connected_pairs);
#    print STDERR "adding edges for queries connected to query $id1; edge count: $edge_count \n" if($edge_count % 100 == 0);
#    for my $id2 (@$conn_to_id1) {
#       $edge_count++;
#       $G->add_edge($id1, $id2);
#    }
#    #  $count++;
# }
# print "# conn pairs: ", scalar keys %connected_pairs, " sum of values:  ", sum(values %connected_pairs), "\n";


# # Get the connected components; each of which is a set of queries which go in the same family.
# my $sum_sizes = 0;
# print STDERR "Now get connected components of graph.\n";
# my %queries_in_biconnected_components = (); # key: core (query) id; value: number of biconn comps it belongs to.
# my @connected_components = $G->connected_components();
# #my @biconnected_components = $G->biconnected_components();
# print STDERR "Done getting connected components. \n";
# print "# Number of connected components: ", scalar @connected_components, "\n";
# #print "# Number of connected components: ", scalar @biconnected_components, "\n";

# for my $ccomponent (@connected_components) { # loop over connected components (i.e. families)
#       $sum_cc_sizes += scalar @$ccomponent;
# }

# # end of loop over connected components.
# #print "# Number of connected components: ", scalar @connected_components, "\n";
# #print "# Number of biconnected components: ", scalar @biconnected_components, "\n";
# print "# vertex count: $count_vertices  edge count: $edge_count  sum of connected component sizes: $sum_sizes \n";

# sub cc_analyze{
#    my $ccomponent = shift;      # (bi) connected component
#    my $query_match = shift;
#    my %id_in_component__q_match_count = (); # keys are ids in family (not just core ids), values are count (of queries that the id matches)
#    my %match_minev = (); # keys are ids in family (not just core), values are min e-values over core ids.
#    my $biggest_q_match_set_size = -1;
#    my $ccsize = scalar @$ccomponent; # @$ccomponent is array of query ids
#    $sum_sizes += $ccsize; # total number of queries which have been dealt with so far.
#    my $sum_id1_matches = 0;
#    for my $core_id (@$ccomponent) { # for each core id ...
#       $queries_in_biconnected_components{$core_id}++;
#       my $q_match_set_size = scalar keys %{$query_match->{$core_id}};
#       $sum_id1_matches += $q_match_set_size;
#       $biggest_q_match_set_size = max($q_match_set_size, $biggest_q_match_set_size);
#       for my $match_id (keys %{$query_match->{$core_id}}) { # increment the count of each id which is a match to it:\
#          my $match_evalue = $query_match->{$core_id}->{$match_id};
#          $id_in_component__q_match_count{$match_id}++;
#          $match_minev{$match_id} = (exists $match_minev{$match_id})? 
#            min($match_minev{$match_id}, $match_evalue) :
#              $match_evalue;
#       }
#    }
#    my @vals = values %id_in_component__q_match_count;
#    my $mult_count = histogram(\@vals);
#    my @sskeys = sort {$b <=> $a} keys %$mult_count;


#    my $fam_size = scalar @vals;
#    my ($overlap90_count, $overlap70_count, $overlap50_count) = (0, 0, 0);
#    my $count_all = 0;
#    for (@sskeys) {
#       #  print "$_: ", $mult_count->{$_}, "; ";
#       if ($_ >= 0.9*$ccsize) { # number of seqs which match at least 90% of the core ids.
#          $overlap90_count += $mult_count->{$_};
#       }
#       if ($_ >= 0.7*$ccsize) { # number of seqs which match at least 80% of the core ids.
#          $overlap70_count += $mult_count->{$_};
#       }
#       if ($_ >= 0.5*$ccsize) { # number of seqs which match at least 50% of the core ids.
#          $overlap50_count += $mult_count->{$_};
#       }
#       $count_all += $mult_count->{$_};
#    }                            # print "\n";
#    my $union_size = scalar keys %match_minev;
#    my $cc_size_all =  scalar keys %id_in_component__q_match_count;
#    printf("%s  ", join(";", sort @$ccomponent)); # the queries of the nodes in connected component 
#    my $cc_edge_count = 0;
#    my $max_cc_edge_count = 0;   # = $ccsize*($ccsize-1)/2;
#    #   for my $cc (@$ccomponent) {
#    #       $cc_edge_count += $G->degree($cc);
#    # #my $pair_id = ($query le $match)? "$query $match" : "$match $query";
#    #    }
#    for my $i (0..scalar @$ccomponent-2) {
#       my $id_i = $ccomponent->[$i];
#       for my $j ($i+1..scalar @$ccomponent-1) {
#          my $id_j = $ccomponent->[$j];
#          my $pair_id = ($id_i le $id_j)? "$id_i $id_j" : "$id_j $id_i";
#          $max_cc_edge_count++;
#          $cc_edge_count++ if(exists $connected_pairs{$pair_id});
#       }
#    }
#    #  $cc_edge_count /= 2;
#    #  my $max_cc_edge_count = $ccsize*($ccsize-1)/2;
#    #print "$cc_edge_count  $max_cc_edge_count \n";
#    printf("%4i %5i  %4.2f %4.2f %4.2f  %4.2f %4.2f %4.2f \n", $ccsize, $cc_size_all, 
#           $sum_id1_matches/($ccsize*$cc_size_all), $biggest_q_match_set_size/$union_size, 
#           ($ccsize > 1)? $cc_edge_count/$max_cc_edge_count : 1,
#           $overlap90_count/$count_all, $overlap70_count/$count_all, $overlap50_count/$count_all);
#    # while(my($k, $v) = each %match_minev){
#    #    print "minev: $k  $v \n";
#    # }
#    for (sort {$match_minev{$a} <=> $match_minev{$b}} keys %match_minev) {

#       printf("    %1s  %8.3g \n", $_ , $match_minev{$_});
#    }
#    print "\n";
# }

# sub ids_to_species{
#    my $ids = shift;          # array ref holding ids
#    my $seqid_species = shift; # keys ids; values corresponding species.
#    my %species = ();

#    for (@$ids) {
#       $species{$seqid_species->{$_}}++;
#    }
#    return \%species;
# }


# sub connected_queries{ # given a query id, return a array ref of other queries connected to it.
#    # but only the ones which weren't previously found and stored in $connected_pairs 
#    my $query = shift;
#    my $q_qAMids = shift;
#    my $id_species = shift;
#    my $min_overlap_fraction = shift;
#    my $connected_pairs = shift;

#    my @connected_to_q = ();
#    for my $qamid ( keys %$q_qAMids) { # %{$q_qAMids->{$query}} ) { # list of ids in qAM subtree of q
#       next if ($query eq $qamid); # don't need to put matches with self into graph.
#       if (exists $q_qAMids->{$qamid}) { # $qamid  is a key of %$q_qAMids, i.e. is itself a query id.
#          my $pair_id = ($query le $qamid)? "$query $qamid" : "$qamid $query";
#          next if (exists $connected_pairs->{$pair_id}); # { # this pair already known to be connected.
#          my $other_query = $qamid;
#          my @set1_ids =  keys  %{$q_qAMids->{$query}};
#          my @set2_ids =  keys  %{$q_qAMids->{$other_query}};
#          my ($size1, $size2) = (scalar @set1_ids, scalar @set2_ids);
#          my $minsize = min($size1, $size2);
#          my $ovrlap = query_pair_connected($query, $other_query, $q_qAMids); # array ref
#          my $ovrlap_size = scalar @$ovrlap;
#          my $set1_species = ids_to_species(\@set1_ids, $id_species);
#          my $set2_species = ids_to_species(\@set2_ids, $id_species);
#          my $overlap_species = ids_to_species($ovrlap, $id_species);
#          my $n_set1_species = scalar keys %$set1_species;
#          my $n_set2_species = scalar keys %$set2_species;
#          my $n_overlap_species = scalar keys %$overlap_species;
#          my $min_species = min($n_set1_species, $n_set2_species);
         
#          if ($ovrlap_size >= 1) {
#             print "$query  $qamid    $size1  $size2  $minsize  $ovrlap_size    $n_set1_species  $n_set2_species  $min_species  $n_overlap_species     ";
#             printf("%8.5f   %8.5f \n", $ovrlap_size/$minsize, $n_overlap_species/$min_species);
#          }
#          #     if (query_pair_connected($query, $other_query, $q_qAMids, $min_overlap_fraction)) {
#          if (  ($ovrlap_size >= $min_overlap_fraction * min($size1, $size2))
#              or
#              ($n_overlap_species >= $min_overlap_fraction * $min_species)  ){
#             push @connected_to_q, $other_query;
#             $connected_pairs->{$pair_id}++;
#          }
#       }
#    }				# end of loop
#    return \@connected_to_q;
# }

# sub query_pair_connected{ # does not use e-value - just uses presence/absence in q_m_href and m_q_href.
#    my $qid1 = shift;
#    my $qid2 = shift;
#    my $q_qAMids = shift;
#    #   my $min_overlap_fraction = shift || 0.9;

#    #   return 0 if(!exists $q_qAMids->{$qid1}->{$qid2}  and  !exists $q_qAMids->{$qid2}->{$qid1}); # require 

#    # my $size1 = scalar keys  %{$q_qAMids->{$qid1}}; # size of set of matches to $qid1
#    # my $size2 = scalar keys  %{$q_qAMids->{$qid2}};
#    # my $Ok = 0;
#    my $ovrlp = overlap($qid1, $qid2, $q_qAMids);
#    # if ($size1 <= $size2) { # overlap defined as fraction >= $min_overlap_fraction
#    #    # of smaller set is in intersection with larger set.
#    #    $Ok = ($ovrlp >= $min_overlap_fraction*$size1);
#    # } else {
#    #    $Ok = ($ovrlp >= $min_overlap_fraction*$size2);
#    # }
#    return $ovrlp;
# }

# # Counts how many are matches to both qid1 and qid2, and how many are matches to qid1, but not qid2.
# sub overlap{ # $qid1 is the id of the query with the smaller set of matches.
#    my $qid1 = shift;
#    my $qid2 = shift;
#    my $q_qAMids = shift;
#    my ($count_both, $count_1not2) = (0, 0);
#    my @intersection_ids = ();
#    for my $qamid (keys %{$q_qAMids->{$qid1}}) {
#       if (exists $q_qAMids->{$qid2}->{$qamid}) { # i.e. $qid2 has $mid1 as a match
#          $count_both++;		# both ev1, ev2 are < threshold
#          push @intersection_ids, $qamid;
#       }
#    }
#    return \@intersection_ids;
# }

# sub histogram{			# histogram an array of integers.
#    my $aref = shift;
#    my %hist = ();
#    for (@$aref) {
#       $hist{$_}++;
#    }
#    return \%hist;
# }

# # sub big_and_union_sizes{
# # my $connected_component = shift;
# # my @cc = @$connected_component; # array of vertices
# # for(@cc){



# sub count_cc_edges{
#    my $conncomp = shift;        # ref to array of vertices
#    my $conn_pairs = shift;	# keys: "id1 id2"
#    my $edge_count = 0;	 # counts the edges in the connected component
#    my @cc = @$conncomp;
#    my $ccsize = scalar @cc;
#    my $max_edges = $ccsize*($ccsize - 1)/2;
#    for (my $i = 0; $i < $ccsize; $i++) {
#       my $id1 = $cc[$i];
#       for (my $j = $i+1; $j < $ccsize; $j++) {
#          my $id2 = $cc[$j];
#          my $idpair = ($id1 le $id2)? "$id1 $id2" : "$id2 $id1";
#          $edge_count++ if(exists $conn_pairs->{$idpair});
#       }
#    }
#    return ($edge_count, $max_edges);
# }



# sub store_gg_info {
#    my $gg_filename   = shift;
#    my %seqid_species = ();
#    my %species_count = ();
#    if ( defined $gg_filename and -f $gg_filename ) {
#       open my $fh_gg, "<", "$gg_filename";
#       while (<$fh_gg>) {
#          my @cols = split( " ", $_ );
#          my $species = shift @cols;
#          $species =~ s/:$//;    # remove final colon if present.
#          $species_count{$species}++;
#          for (@cols) {
#             my $id = $_;
#             $id =~ s/[|]/_/g;   # change pipes to underscores.
#             $seqid_species{$id} = $species;
#          }
#       }
#    }                # done storing gg_file info in hash %seqid_species
#    return (\%seqid_species, \%species_count);
# }
