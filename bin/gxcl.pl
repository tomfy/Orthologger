#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw( min max sum );
no lib '/home/tomfy/bin';
no lib '/home/tomfy/Orthologger/bin';
no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}

use lib $libdir;

use Getopt::Long;
use CXGN::Phylo::BasicTree;
use Hash::Ordered;
use TomfyMisc qw ( increment_hash  add_hashes  format_newick_species_info median  store_gg_info );
use Grouper;

=pod
uta: unrooted tree analyzer
analyze a set of trees in an unrooted fashion.
the general idea is that we are interested in finding subtrees
each of which contains at least one sequence belonging to a set of
species, and no sequences belonging to
some other set of species.
e.g. we could require at least one sequence from a coffee species
(say C. arabica, C. canephora or C. eugenioides), and zero sequences of
Populus or Arabidopsis.
For each coffee sequence, root the tree near that sequence
then work down through the tree, until reaching a node where both
subtrees contain disallowed species. The branch above that node defines
a bipartition of the tree, with one part the maximal subtree containing the
starting leaf, and no disallowed species.

-in  <pattern specifying a set of files each containing a single newick expression>
-gr  <filename of file specifying a set of species, all leaves of which will be investigated, and specifying a set of disallowed species.>

Usage example:
uta.pl -in '*.newick' -gr pgroup  -min_max 6  -max_to_med 4  >  all_6_4.out

=cut

my $D2group = 'brassicaceae,3';
my $Ogroupname = 'to_output';
{
   # Defaults:
   my $input_pattern = undef;
   my $gg_filename   = undef;
   my $groups = 'groups'; # this is a default filename specifying the groups of  species.
   my $min_sgroup_sequences = 1;
   my $min_max_sgroup_expansion = 1;
   my $max_to_median_factor = 1;
   my $n_disallowed_limit = 0;
   my $analysis_type = 1;
   my $S1group = 'dicots,6,2';
   my $S2group = 'monocots,4,1';
   my $D1group = 'basals,1';    # disallowed groupname, max allowed
   #   my $Fgroupname = 'nonams,1';
   

   GetOptions(
              'input_pattern=s' => \$input_pattern, # e.g. '*.newick'
              'groups=s' => \$groups, # filename or string (e.g. 'Coffea_arabica,Coffea_canephora') specifying positive species.
              'dis_limit=i' => \$n_disallowed_limit, # maximum number of disallowed sequences in subtree (either 0 or 1)
              'min_sequences=i' => \$min_sgroup_sequences, # min total sgroup sequences in no-disallowed subtree.
              'min_max_exp=i' => \$min_max_sgroup_expansion, # require at least one sgroup sp has at least this many sequences in subtree.
              'max_to_median=f' => \$max_to_median_factor, # require that max number of sgroup sequences is at least this many times median.
              'ggfile=s'    => \$gg_filename, # defines species-sequence association (not needed if species info present in newick.)
              'analysis_type=s' => \$analysis_type,
              's1=s' => \$S1group,
              's2=s' => \$S2group,
              'd1=s' => \$D1group,
             );

   ######   get the input filenames (should be newick format)  ######
   die "input pattern is undefined.\n" if(!defined $input_pattern);
   my $input_filenames_str = `ls $input_pattern`;
   die "No files found matching pattern: $input_pattern \n" if($input_filenames_str eq '');
   my @input_filenames = split(" ", $input_filenames_str);

   ######  store the gene-genome information  ######
   my ($gg_seqid_species, $gg_species_count) = (defined $gg_filename)? store_gg_info($gg_filename) : (undef, {});
   my @sspecies = sort {$a cmp $b}  keys %$gg_species_count; # sorted species

   my $grouper = Grouper->new($groups);
   print $grouper->get_group_species_string();

   print STDERR "# Will analyze ", scalar @input_filenames, " files.\n";
   for my $input_filename (@input_filenames) {
      print STDERR "$input_filename\n";
      open my $fh_in, "<", "$input_filename" or die "Couldn't open $input_filename for reading.\n";
      my $input_newick = '';
      while (<$fh_in>) {
         $input_newick .= $_;
      }
      $input_newick =~ s/\s//g;
      $input_newick =~ s/;?\s*$//; # remove final ; and whitespace if present

      if ($input_newick) {
         my ($newick, $seqid_species) =  format_newick_species_info($input_newick, 1);
         my $the_tree = make_tree($newick);

         if ($analysis_type  eq  1) { # max subtree without disalloweds (or with <= 1).
            my $output_strings = analyze_tree($the_tree, $grouper, undef, $seqid_species,
                                              $n_disallowed_limit, $min_sgroup_sequences, $min_max_sgroup_expansion, $max_to_median_factor);
            for my $outstr (@$output_strings) { # loop over 'good' subtrees, outputting info for each
               print "$input_filename  $outstr \n";
            }
         } elsif ($analysis_type eq 2) {
            my @output_strings = analyze_tree_2($the_tree, $grouper, $seqid_species, $S1group, $S2group, $D1group); #, $Fgroupname);
            for my $outstr (@output_strings) {
               print "$input_filename  $outstr\n";
            }
         } else {               # whole tree
            my $sgroup_name = $grouper->get_ith_groupname(0);
            $grouper->group_species_counts($the_tree->get_root()->get_implicit_species()); #, $seqid_species);
            my ($nsppres, $med_nseq, $max_nseq, $total_nseq) = $grouper->get_group_summary_info($sgroup_name);
            if ( ($nsppres >= $min_sgroup_sequences) and ($max_nseq >= $min_max_sgroup_expansion)  and  ($max_nseq >= $max_to_median_factor * $med_nseq) ) {
               print "$input_filename    ";
               #           print STDERR "group names:  ", join("  ", $grouper->get_groupnames()), "\n";
               for my $grpname ($grouper->get_groupnames()) {
                  next if($grpname eq 'disallowed');
                  printf("%1i ", $grouper->get_n_species_in_group($grpname) );
                  print join("  ", $grouper->group_speciescount_strings($grpname)), "     ";
               }
               my %seqid_presence = (map(($_ => 1), @{$the_tree->get_root()->get_implicit_names()}));
               print $grouper->ids_present_string($seqid_species, \%seqid_presence, $sgroup_name),  "\n";
            }
         }
      } else {                  # the newick string is empty
         print STDERR "File $input_filename does not contain a newick string.\n";
      }
   }                            # end of loop over input files
   print STDERR "# Done.\n";
}
# end of main


##############  Subroutines  ################
sub get_species_group_leaves{
   my $tree = shift;
   my $group = shift;    # Hash::Ordered of species (values are all 1)
   my @selected_leaves = ();
   for my $a_leaf ($tree->get_leaves()) {
      my $leaf_species = $a_leaf->get_species();
      #    print "leaf name, species:  ", $a_leaf->get_name(), "  ", $leaf_species, "\n";
      #  if (exists $group->{$leaf_species}) {
      if (defined $group->get($leaf_species)) { # if leaf species is in group, add leaf node to @selected_leaves
         #    print "$leaf_species is in group.\n";
         push @selected_leaves, $a_leaf;
      }
   }
   return \@selected_leaves;
}

sub make_tree{
   my $the_input_newick = shift;
   my $parser = CXGN::Phylo::Parse_newick->new( $the_input_newick, 0 );
   my $tree = $parser->parse( CXGN::Phylo::BasicTree->new() );
   $tree->make_binary();
   $tree->impose_branch_length_minimum(); # using default min branch length
   $tree->show_newick_attribute("species");
   $tree->set_show_standard_species(0);
   $tree->get_root()->recursive_implicit_names();
   $tree->get_root()->recursive_implicit_species();
   return $tree;
}

sub analyze_tree{
   my $tree = shift;
   my $grouper = shift;           # Grouper object
   my $select_group_name = shift; # if undef, use 0th group of grouper obj.
   my $sequenceid_species = shift;

   my $n_disallowed_limit = shift;
   my $min_sgroup_seqs = shift;
   my $min_max_sgroup_exp = shift;
   my $sgroup_max_to_median = shift;
   my $disallowed_limit = shift // 0;

   my $done_ids = {};        # keys are ids which are already analyzed
   my @output_strings = ();

   $select_group_name //= $grouper->get_ith_groupname(0); # specifies group upon which selection will be based. Default is first group.
   my $leaves_to_analyze = get_species_group_leaves( $tree, $grouper->get_group($select_group_name) );

   my $count_leaves_done = 0;
   my $count_leaves_analyzed = 0;
   my $reroot_count = 0;
   for my $the_leaf (@$leaves_to_analyze) { # loop over leaves of interest
      $count_leaves_done++;
     
      my $sequence_id = $the_leaf->get_name();
      #   print STDERR "  $sequence_id \n";
      next if($done_ids->{$sequence_id}); # the subtree containing this one has already been found - skip.
      $reroot_count++;
      #  print STDERR "     rerooting near $sequence_id. $reroot_count \n";
      ######  Reroot near the leaf  ######
      $count_leaves_analyzed++;
      my $seqid_presence = {$sequence_id => 1}; # this hash holds the ids in the maximal part of tree containing query and only pgroup sequences.

      $tree->reset_root_to_point_on_branch($the_leaf, 0.5*$the_leaf->get_branch_length());
      $tree->get_root()->recursive_implicit_names(); #

      my @children = $tree->get_root()->get_children();
      die "Root node has ", scalar @children, " children; should have exactly 2.\n" if(scalar @children != 2);
      my $next_node = $children[1];
      my $the_subtree_string = "[" . $the_leaf->get_name() . " ";
      my $subtree_string_other_part = '';
      my $n_dis_so_far = 0;
     
      while (1) {
         @children = $next_node->get_children();
         my $n_children = scalar @children;
         #$Lspecies_count : keys: species, values: number of leaves in subtree of that species.
         if ($n_children == 2) {
            my ($Lchild, $Rchild) = @children;
            my ($Lspecies_count, $Lcat_spcount, $Lcat_leafcount, $Lid_presence) =  $grouper->species_and_group_counts($Lchild->get_implicit_names(), $sequenceid_species);
            my ($Rspecies_count, $Rcat_spcount, $Rcat_leafcount, $Rid_presence) =  $grouper->species_and_group_counts($Rchild->get_implicit_names(), $sequenceid_species);
            my $L_n_dis = $Lcat_leafcount->{'disallowed'} // 0; # number of leaves of disallowed species found in L subtree
            my $R_n_dis = $Rcat_leafcount->{'disallowed'} // 0; # number of leaves of disallowed species found in R subtree
            my $L_n_leaves = scalar keys %$Lid_presence;
            my $R_n_leaves = scalar keys %$Rid_presence;
            # max subtree with <= $n_disallowed_limit disallowed leaves
            if ($L_n_dis + $R_n_dis + $n_dis_so_far <= (($n_disallowed_limit == 1)? 2 : 0)) { # add both, and done
               $seqid_presence = increment_hash($seqid_presence, $Lid_presence);
               $seqid_presence = increment_hash($seqid_presence, $Rid_presence);
               $subtree_string_other_part = $Lchild->get_implicit_names()->[0] . " ()";
               last;
            } elsif ($L_n_dis + $n_dis_so_far <= $n_disallowed_limit) { # add L, continue
               $seqid_presence = increment_hash($seqid_presence, $Lid_presence);
               $next_node = $Rchild;
               $subtree_string_other_part = $Lchild->get_implicit_names()->[0];
               $n_dis_so_far += $L_n_dis;
            } elsif ($R_n_dis + $n_dis_so_far <= $n_disallowed_limit) { # add R, continue
               $seqid_presence = increment_hash($seqid_presence, $Rid_presence);
               $next_node = $Lchild;
               $subtree_string_other_part = $Rchild->get_implicit_names()->[0];
               $n_dis_so_far += $R_n_dis;
            } else {            # add neither, done
               $subtree_string_other_part .= " (" . $next_node->get_implicit_names()->[0]. ")";
               last;
            }
         } elsif ($n_children == 0) { # to get here next_node must be a single nonAM leaf -- done
            last;
         } elsif ($n_children > 2) {
            die "node: ", $next_node->get_name(), " has ", scalar @children, " children (not binary tree).\n"
         } else {
            die "node: ", $next_node->get_name(), " has ", scalar @children, " children. ???.\n"
         }
      }                         # end of while(1) loop

      my $species_count = ();
      for my $k (keys %$seqid_presence) {
         $done_ids->{$k} = 1;
         $species_count->{$sequenceid_species->{$k}}++;
      }

      $grouper->group_species_counts($species_count);
      my ($nsppres, $med_nseq, $max_nseq, $total_nseq) = $grouper->get_group_summary_info( $select_group_name );
      if ( ($total_nseq >= $min_sgroup_seqs) and
           ($max_nseq >= $min_max_sgroup_exp) and
           ($max_nseq >= $sgroup_max_to_median*$med_nseq)
         ) {
         my $outstring = '';
         for my $a_grpname ($grouper->get_groupnames()) {
            next if($a_grpname eq 'disallowed');
            $outstring .= sprintf("%1i  ", scalar $grouper->get_group($a_grpname)->keys());
            $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($a_grpname)) );
         }
         $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $select_group_name);
         push @output_strings, $outstring;
      }
   }                            # loop over leaves to analyze
   print STDERR "     rerooted $reroot_count times.\n"; #exit;
   $tree->impose_branch_length_minimum(1);
   $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
   return \@output_strings;
}

sub analyze_tree_2{
   my $tree = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $S1group = shift;
   my $S2group = shift;
   my $D1group = shift;
   #  my $Fgroupname = shift;
   my $node = $tree->get_root();
   return recursive_is_node_speciation($node, $grouper, $sequenceid_species, $S1group, $S2group, $D1group); #, $Fgroupname);
}

sub recursive_is_node_speciation{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $S1group = shift;
   my $S2group = shift;
   my $D1group = shift;
   my ($is_spec, $which, $outstring) = is_node_speciation($node, $grouper, $sequenceid_species, $S1group, $S2group, $D1group);
   my @strings = ($is_spec)? ($outstring) : ();
   if (1==1) {
      if (1 or (! ($is_spec and ($which eq 'LR')))) { # if child subtrees (i.e. LR) look like speciation don't descend to child nodes.
         my @children = $node->get_children();
         for (@children) {
            push @strings, recursive_is_node_speciation($_, $grouper, $sequenceid_species, $S1group, $S2group, $D1group);
         }
      }
   }
   return @strings;
}

sub two_subtrees_speciation{ 
   # given $G1_spcount, and $G2_spcount, the counts of species in each group
   # present in subtrees 1 and 2;
   # check whether it looks like a speciation, with one subtree having many 
   # species from Lgroup, and few from Rgroup and vice versa
   # and both subtrees have few from Bgroup
   my $G1_spcount = shift;
   my $G2_spcount = shift;
   #  my $Lgroupname = shift;
   my $S1group = shift;
   my $S2group = shift;
   my $D1group = shift;
   my $which2 = shift;
   #   print "L,R groups:  $Lgroup  $Rgroup \n";
   my ($Lgroupname, $Lgrp_min, $Lgrp_max) = split(/[,\s]+/, $S1group); # (6, 2) the min in one subtree, and the max in the other
   my ($Rgroupname, $Rgrp_min, $Rgrp_max) = split(/[,\s]+/, $S2group); #(4, 1);
   my ($Bgroupname, $Bgrp_max) = split(/[,\s]+/, $D1group); #1; # both subtrees should have at most this many 'basal' species.
   my ($Dgroupname, $Dgrp_max) = split(/[,\s]+/, $D2group); #1; # both subtrees should have at most this many disallowed species.

   my $G1_B = $G1_spcount->{$Bgroupname} // 0;
   my $G2_B = $G2_spcount->{$Bgroupname} // 0;
   my $G1_D = $G1_spcount->{$Dgroupname} // 0;
   my $G2_D = $G2_spcount->{$Dgroupname} // 0;
   #   print "G1_B, G2_B: $G1_B $G2_B \n";
   if (                         # both G1 and G2 have few Bgroup
       ($G1_B <= $Bgrp_max)
       and
       ($G2_B <= $Bgrp_max)
       and
       ($G1_D <= $Dgrp_max)
       and
       ($G2_D <= $Dgrp_max)
      ) {
      my $G1_L = $G1_spcount->{$Lgroupname} // 0;
      my $G1_R = $G1_spcount->{$Rgroupname} // 0;
      my $G2_L = $G2_spcount->{$Lgroupname} // 0;
      my $G2_R = $G2_spcount->{$Rgroupname} // 0;
      #    print "  $which2    $G1_L  $G1_R  $G1_B   $G2_L  $G2_R  $G2_B \n"; 
      if ( (
            # G1 has many Lgroup and few Rgroup
            ($G1_L >= $Lgrp_min)
            and
            ($G1_R <= $Rgrp_max)
            and                 # G2 has many Rgroup and few Lgroup, 
            ($G2_L <= $Lgrp_max)
            and
            ($G2_R >= $Rgrp_min)
           )
           or (
               # G2 has many Lgroup and few Rgroup
               ($G2_L >= $Lgrp_min)
               and
               ($G2_R <= $Rgrp_max)
               and              # G1 has many Rgroup and few Lgroup, 
               ($G1_L <= $Lgrp_max)
               and
               ($G1_R >= $Rgrp_min)
              )
         ) {
         return (1 == 1); # true (looks like a speciation (approximately))
      }
   } else {
      return (1 == 0);        # false (doesn't look like a speciation)
   }
}


sub is_node_speciation{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $S1group = shift;
   my $S2group = shift;
   my $D1group = shift;
   #  my $group = shift;

   my $which = undef;
   my ($Lspecies_count, $Lcat_spcount, $Lid_presence, 
       $Rspecies_count, $Rcat_spcount, $Rid_presence,
       $Aspecies_count, $Acat_spcount, $Aid_presence) 
     = $grouper->three_subtrees_species_counts($node, $sequenceid_species);
   my $Lgroupname = ($S1group =~ /\s*([^,\s]+)/)? $1 : '';
   my $Rgroupname = ($S2group =~ /\s*([^,\s]+)/)? $1 : '';
   my $Bgroupname = ($D1group =~ /\s*([^,\s]+)/)? $1 : '';
   my $species_count;
   my $outstring = '';
my $seqid_presence;

   if (two_subtrees_speciation($Lcat_spcount, $Rcat_spcount, $S1group, $S2group, $D1group, "LR")) {
      print STDERR "L R  ", 
        $Lcat_spcount->{$Lgroupname} // 0, "  ", $Lcat_spcount->{$Rgroupname} // 0, "  ", $Lcat_spcount->{$Bgroupname} // 0, "      ",
          $Rcat_spcount->{$Lgroupname} // 0, "  ", $Rcat_spcount->{$Rgroupname} // 0, "  ", $Rcat_spcount->{$Bgroupname} // 0, "\n";
      $which = 'LR';
      $species_count = add_hashes($Lspecies_count, $Rspecies_count);
      $seqid_presence = add_hashes($Lid_presence, $Rid_presence);
      $grouper->group_species_counts($species_count);
   } elsif (two_subtrees_speciation($Rcat_spcount, $Acat_spcount, $S1group, $S2group, $D1group, "RA")) {
      print STDERR "R B  ", 
        $Rcat_spcount->{$Lgroupname} // 0, "  ", $Rcat_spcount->{$Rgroupname} // 0, "  ", $Rcat_spcount->{$Bgroupname} // 0, "      ",
          $Acat_spcount->{$Lgroupname} // 0, "  ", $Acat_spcount->{$Rgroupname} // 0, "  ", $Acat_spcount->{$Bgroupname} // 0, "\n";
      $which = 'RB';
      $species_count = add_hashes($Rspecies_count, $Aspecies_count);
$seqid_presence = add_hashes($Rid_presence, $Aid_presence);
      $grouper->group_species_counts($species_count);
   } elsif (two_subtrees_speciation($Acat_spcount, $Lcat_spcount, $S1group, $S2group, $D1group, "AL")) {
      print STDERR "B L  ", 
        $Acat_spcount->{$Lgroupname} // 0, "  ", $Acat_spcount->{$Rgroupname} // 0, "  ", $Acat_spcount->{$Bgroupname} // 0, "      ",
          $Lcat_spcount->{$Lgroupname} // 0, "  ", $Lcat_spcount->{$Rgroupname} // 0, "  ", $Lcat_spcount->{$Bgroupname} // 0, "\n";
      $which = 'BL';
      $species_count = add_hashes($Aspecies_count, $Lspecies_count);
      $seqid_presence = add_hashes($Aid_presence, $Lid_presence);
 #     while(my($k,$v) = each %$species_count){ print "$k:$v  "; } print "\n";
      $grouper->group_species_counts($species_count);
   } else {                     # doesn't look like a speciation
      return (1 == 0, $which, $outstring);
   }
   for my $gname ($Lgroupname, $Rgroupname) {
      $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
      $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
   }
   $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);
   return (1 == 1, $which, $outstring);
}


sub LR_ids{   # get one id from L child subtree, and one for R subtree
   # and return string "[$Lid $Rid]"
   # or (if node is leaf): "[leaf:  Leafid]"
   my $node = shift;
   my $subtree_str = '[';
   if ($node->is_leaf()) {
      $subtree_str = "[leaf " . $node->get_name() . "]";
   } else {
      my @children = $node->get_children();
      for my $i (0..1) {   # loop over children 0 and 1 (i.e. L and R)
         my @implicit_names = @{$children[$i]->get_implicit_names()};
         $subtree_str .= $implicit_names[0] . " ";
      }
      $subtree_str =~ s/\s+$/]/;
      #$subtree_str .= "]";
   }
   return $subtree_str;
}

#######################################
# unused
#######################################

sub max_species_leafcount_of_category{
   my $species_leafcount = shift;
   my $species_category = shift;
   my $category = shift;        # e.g. 'disallowed'
   my $max_species_leafcount = 0;
   while (my($sp, $count) = each %$species_leafcount) {
      if ($species_category->{$sp} eq $category
          and
          $count > $max_species_leafcount) {
         $max_species_leafcount = $count;
      }
   }
   return $max_species_leafcount;
}

sub species_category_id_info_string{
   my $ids_in_subtree = shift ;
   my $seqid_sp = shift;
   my $group = shift;         # Hash::Ordered; keys are species names.
   my $output_the_seqids = shift || 0;
   my $output_string = '';
   my %groupspecies_counts = ();
   my @group_ids = ();
   for my $an_id (@$ids_in_subtree) {
      my $sp = $seqid_sp->{$an_id};
      if (defined $group->get($sp)) {
         $groupspecies_counts{$sp}++;
         push @group_ids, $an_id
      }
   }
   $output_string .= sprintf("%1i %1i  ", scalar keys %groupspecies_counts, scalar @group_ids); # number of species in group present in subtree
   my @sorted_psp_counts = map($groupspecies_counts{$_} // 0, $group->keys());
   for my $count (@sorted_psp_counts) {
      $output_string .= sprintf("%1i ", $count);
   }
   my $seqid_string = sprintf(" [%s]", join(";", sort @group_ids)) if($output_the_seqids);
   return ($output_string, $seqid_string);
}

sub species_counts{
   my $node = shift;

   my $sp_count = {};
   my $species = $node->get_implicit_species();
   for my $sp (@$species) {
      $sp_count->{$sp}++;
   }
   return $sp_count;
}
