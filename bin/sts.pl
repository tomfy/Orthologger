#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw( min max sum );

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}

use lib $libdir;

use Getopt::Long;
use CXGN::Phylo::BasicTree;
use Hash::Ordered;
use TomfyMisc qw ( increment_hash  add_hashes  format_newick_species_info median  store_gg_info );
use Grouper;

my $LROK_descend_into_children = 0;

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

=cut

my $min_s_species_present = 12;
my $bigger_than_any_tree = 1e100;

{
   my $cmmnd =  qx/ps -o args $$/;
   $cmmnd =~ s/\n/  /g;
   print "#  $cmmnd \n";

   # Defaults:
   my $input_pattern = undef;
   my $gg_filename   = undef;
   my $groups_file = 'groups'; # this is a default filename specifying the groups of  species.

   my $analysis_type = 1;

   my $Sp_grp1 = 'dicots,6,2';
   my $Sp_grp2 = 'monocots,4,1';
   my $split = 'dicots,6,2;monocots,4,1'; # for use by an2, same info as Sp_grp1, $Sp_grp2.

   my $S12sum_range = '13,10000';
   my $D1group = 'basals,1';    # disallowed groupname, max allowed
   my $D2group = 'brassicaceae,30';
my $Dgroups = undef; # e.g. 'basals;non_ams'  - groups which are strictly not disallowed.
   my $D12sum_range = '0,2';

   my $include_spname_in_seqid = 0;
 #  my $S3group = 'monocots';

   my $speciescount_ranges_str = ''; # e.g. 'erysimum,15,20;brassicaceae,0,5' from 15 to 20 erysimum sp. present, max of 5 brass7 sp. present.
   my $leafcount_ranges_str = ''; # e.g. 'erysimum,30,10000;brassicaceae7,0,9' as many groups as you want, with min and max numbers of sequences of each specified


   GetOptions(
              'input_pattern=s' => \$input_pattern, # e.g. '*.newick'
              'groups_file=s' => \$groups_file, # name of file defining groups of species.

              'ggfile=s'    => \$gg_filename, # defines species-sequence association (not needed if species info present in newick.)

              'analysis_type=s' => \$analysis_type,

              's1=s' => \$Sp_grp1,
              's2=s' => \$Sp_grp2,
              'split=s' => \$split,

              's12sum_range=s' => \$S12sum_range, # e.g. '12,1000'
#              's3=s' => \$S3group,
              'd1=s' => \$D1group,
              'd2=s' => \$D2group,
              'd12sum_range=s' => \$D12sum_range,
              'disallowed_groups=s' => \$Dgroups, # e.g. 'basals,non_ams'

              'speciescount_ranges=s' => \$speciescount_ranges_str, 
              'leafcount_ranges=s' => \$leafcount_ranges_str,

              'prepend_sp=i' => \$include_spname_in_seqid,
             );

   print "xxx D1group $D1group   D2group $D2group \n";
   print "S12sum_range: $S12sum_range.\n";

   ($Sp_grp1, $Sp_grp2) = split(";", $split);
   print "Sp_grp1: $Sp_grp1, Sp_grp2: $Sp_grp2 \n";

   my $group_speciescount_ranges = group_ranges($speciescount_ranges_str);
   my $group_leafcount_ranges = group_ranges($leafcount_ranges_str);

   print STDERR 'D1group: ', $D1group, "\n";

$Dgroups =~ s/ //g;
my @disallowed_groups = split(/[,;]/, $Dgroups);

   ######   get the input filenames (should be newick format)  ######
   die "input pattern is undefined.\n" if(!defined $input_pattern);
   my $input_filenames_str = `ls $input_pattern`;
   die "No files found matching pattern: $input_pattern \n" if($input_filenames_str eq '');
   my @input_filenames = split(" ", $input_filenames_str);

   ######  store the gene-genome information  ######
   my ($gg_seqid_species, $gg_species_count) = (defined $gg_filename)? store_gg_info($gg_filename) : (undef, {});
   my @sspecies = sort {$a cmp $b}  keys %$gg_species_count; # sorted species

   my $grouper = Grouper->new($groups_file);
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

         my ($newick, $seqid_species) =  format_newick_species_info($input_newick, 1, $include_spname_in_seqid);
         my $the_tree = make_tree($newick);

         if ($analysis_type  eq  1) { # max subtree without disalloweds (or with <= 1).
            die "analysis type 1 not implemented.\n";
            # my $output_strings = analyze_tree_1($the_tree, $grouper, undef, $seqid_species,
            #                                     $n_disallowed_limit, $min_sgroup_sequences, $min_max_sgroup_expansion, $max_to_median_factor);
            # for my $outstr (@$output_strings) { # loop over 'good' subtrees, outputting info for each
            #    print "$input_filename  $outstr \n";
            # }
         } elsif ($analysis_type eq 2) {
                 print STDERR "analysis type 2, before analysis: \n";
                 print "CCC: $Sp_grp1 $Sp_grp2 $D1group $D2group  $S12sum_range  $D12sum_range \n";
            my @which_output_strings = analyze_tree_2($the_tree, $grouper, $seqid_species, 
                                                      $Sp_grp1, $Sp_grp2, $D1group, $D2group, 
                                                      $S12sum_range, $D12sum_range);
            #  print "done analyzing tree.\n";
            for my $which_outstring (@which_output_strings) {
               for my $wh ('LR', 'RB', 'BL') {
                  if (exists $which_outstring->{$wh}) {
                     my $outstr = $which_outstring->{$wh} // '';
                     print "$input_filename    $outstr"; # \n";
                  }
               }
            }
 #        } elsif ($analysis_type eq 3) { # just look for subtrees with as many from Sp_grp1 as possible, with none from D1group.
 #           analyze_tree_3($the_tree, $grouper, $seqid_species, $Sp_grp1, $D1group, $input_filename);
         } elsif ($analysis_type eq 4) { # look for loss of one group after duplication, e.g. erysimum present in both subtrees, other brass. in only one.
            #print STDERR "analysis type 4: \n";
            analyze_tree_4($the_tree, $grouper, $seqid_species,
                           $group_speciescount_ranges,
                           $group_leafcount_ranges,
                           \@disallowed_groups,
                           $input_filename);
 #        } elsif ($analysis_type eq 5) {
 #           analyze_tree_5($the_tree, $grouper, $seqid_species, $Sp_grp1, $Sp_grp2, $D1group, $input_filename);
         } else {               # whole tree
            # my $sgroup_name = $grouper->get_ith_groupname(0);
            # $grouper->group_species_counts($the_tree->get_root()->get_implicit_species()); #, $seqid_species);
            # my ($nsppres, $med_nseq, $max_nseq, $total_nseq) = $grouper->get_group_summary_info($sgroup_name);
            # if ( ($nsppres >= $min_sgroup_sequences) and ($max_nseq >= $min_max_sgroup_expansion)  and  ($max_nseq >= $max_to_median_factor * $med_nseq) ) {
            #    print "$input_filename    ";
            #    #           print STDERR "group names:  ", join("  ", $grouper->get_groupnames()), "\n";
            #    for my $grpname ($grouper->get_groupnames()) {
            #       next if($grpname eq 'disallowed');
            #       printf("%1i ", $grouper->get_n_species_in_group($grpname) );
            #       print join("  ", $grouper->group_speciescount_strings($grpname)), "     ";
            #    }
            #    my %seqid_presence = (map(($_ => 1), @{$the_tree->get_root()->get_implicit_names()}));
            #    print $grouper->ids_present_string($seqid_species, \%seqid_presence, $sgroup_name),  "\n";
            # }
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
   my $grouper = shift;
   my $groupnames = shift; # array ref

   my @selected_leaves = ();
   for my $a_leaf ($tree->get_leaves()) {
      my $leaf_species = $a_leaf->get_species();
      my $leaf_is_disallowed = 0;
      for my $group_name (@$groupnames) { # find if leaf belongs to any of the disallowed groups
      #   print STDERR "DDD: group name: $group_name  leaf name, species:  ", $a_leaf->get_name(), "  ", $leaf_species, "\n";
         $leaf_is_disallowed = 1 if ($grouper->get_groupname_species()->get($group_name)->exists($leaf_species));
      }
      push @selected_leaves, $a_leaf if($leaf_is_disallowed);
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



sub analyze_tree_2{ # look for speciation-like nodes: i.e. of three subtrees one is mostly Sp_grp1, another mostly Sp_grp2,
   # and the 3rd can have anything. Small admixtures of wrong species (e.g. S1 in S2 subtree, or D1 or D2) are allowed.
   my $tree = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $Sp_grp2 = shift;
   my $D1group = shift;
   my $D2group = shift;
   my $S12sum_range = shift;
   my $D12sum_range = shift;
   my $node = $tree->get_root();
print "In an.._2.\n";
   print "D12 sum range: $D12sum_range \n";
print STDERR "BBBB D1group: $D1group \n"; # exit;
   my @return_array = recursive_is_node_speciation($node, $grouper, $sequenceid_species, $Sp_grp1, $Sp_grp2, $D1group, $D2group, 
                                                   $S12sum_range, $D12sum_range); #, $Fgroupname);
   $tree->impose_branch_length_minimum(1);
   $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
   return @return_array;
}

sub recursive_is_node_speciation{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $Sp_grp2 = shift;
   my $D1group = shift;
   my $D2group = shift;
   my $S12sum_range = shift;
   my $D12sum_range = shift;
   print "In rec...node_specieation. before is_node_spec..\n";
   my $which_outstring =  is_node_speciation($node, $grouper, $sequenceid_species, 
                            $Sp_grp1, $Sp_grp2, $D1group, $D2group,
                                             $S12sum_range, $D12sum_range);

   my $LROK = exists $which_outstring->{'LR'};
   my @whichoutstrings = (keys %$which_outstring)? ($which_outstring) : ();

   if ($LROK_descend_into_children  or  !$LROK) {
      #    if (1 or (! ($is_spec and ($which eq 'LR')))) { # if child subtrees (i.e. LR) look like speciation don't descend to child nodes.
      my @children = $node->get_children();
      for (@children) {
         push @whichoutstrings, recursive_is_node_speciation($_, $grouper, $sequenceid_species, 
                                                             $Sp_grp1, $Sp_grp2, $D1group, $D2group, 
                                                             $S12sum_range, $D12sum_range); #, $Fgroupname);
      }
      #    }
   }
   return @whichoutstrings;
}

sub two_subtrees_speciation{ 
   # given $G1_spcount, and $G2_spcount, the counts of species in each group
   # present in subtrees 1 and 2;
   # check whether it looks like a speciation, with one subtree having many 
   # species from Lgroup, and few from Rgroup and vice versa
   # and both subtrees have few from Bgroup
   my $G1_spcount = shift;     # the species and counts in one subtree
   my $G2_spcount = shift; # the species and counts in the other subtree
   #  my $Lgroupname = shift;
   my $Sp_grp1 = shift; # e.g. 'dicots,8,2' a group, and the min number of species of this group in one subtree, the max number of species in the other.
   my $Sp_grp2 = shift; # e.g. 'monocots,5,1'
   my $D1group = shift; # e.g. 'basals,1' a group, and the max number of species of that group to allow in each of the two subtrees being tested.
   my $D2group = shift;
   my $which2 = shift;

   #print STDERR "groups: ", join(" ", keys %$G1_spcount), "\n";
print "Sp_grp1: $Sp_grp1   Sp_grp2: $Sp_grp2 \n";
   print "D1group: $D1group   D2group: $D2group \n";

   #   print "L,R groups:  $Lgroup  $Rgroup \n";
   my ($Lgroupname, $Lgrp_min, $Lgrp_max) = split(/[,\s]+/, $Sp_grp1); # (6, 2) the min in one subtree, and the max in the other
   my ($Rgroupname, $Rgrp_min, $Rgrp_max) = split(/[,\s]+/, $Sp_grp2); #(4, 1);
   my ($Bgroupname, $Bgrp_max) = split(/[,\s]+/, $D1group); #1; # both subtrees should have at most this many 'basal' species.
   my ($Dgroupname, $Dgrp_max) = split(/[,\s]+/, $D2group); #1; # both subtrees should have at most this many disallowed species.
   print STDERR "AAAA: $D1group  $Bgroupname  [$Bgrp_max] \n";

   my $G1_B = $G1_spcount->{$Bgroupname} // 0; # number of 'basal' disallowed species in subtree 1
   my $G2_B = $G2_spcount->{$Bgroupname} // 0; # number of 'basal' disallowed species in subtree 2
   my $G1_D = $G1_spcount->{$Dgroupname} // 0; # number of other disallowed species in subtree 1
   my $G2_D = $G2_spcount->{$Dgroupname} // 0; # number of other disallowed species in subtree 2
   my $G1_L = '-';
   my $G2_L = '-';
   my $G1_R = '-';
   my $G2_R = '-';
   #   print STDERR "$which2   G1_B, G2_B: $G1_B $G2_B  G1_D, G2_D:  $G1_D $G2_D   $Bgrp_max $Dgrp_max \n";
   if (                         # both G1 and G2 have few Bgroup
       ($G1_B <= $Bgrp_max)
       and
       ($G2_B <= $Bgrp_max)
       and
       ($G1_D <= $Dgrp_max)
       and
       ($G2_D <= $Dgrp_max)
      ) {
      $G1_L = $G1_spcount->{$Lgroupname} // 0; # number of Lgrp species in subtree 1
      $G1_R = $G1_spcount->{$Rgroupname} // 0; # number of Rgrp species in subtree 1
      $G2_L = $G2_spcount->{$Lgroupname} // 0; # number of Lgrp species in subtree 2
      $G2_R = $G2_spcount->{$Rgroupname} // 0; # number of Rgrp species in subtree 2
      print "ZZZ  $which2    $G1_L  $G1_R  $G1_B   $G2_L  $G2_R  $G2_B \n"; 
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
             print  "OK: $which2   $G1_L $G1_R   $G2_R $G2_L   $G1_B $G2_B   $G1_D $G2_D \n" if(1 or $which2 =~ /B/); 
         return (1 == 1); # true (looks like a speciation (approximately))
      } else {
             print  "NG2: $which2   $G1_L $G1_R   $G2_R $G2_L   $G1_B $G2_B   $G1_D $G2_D \n" if(1 or $which2 =~ /B/); 
         #   print STDERR "NG: $which2   $G1_L $G2_L   $G1_R $G2_R   $G1_B $G2_B   $G1_D $G2_D \n" if($which2 =~ /A/); 
         return (1 == 0);
      }
   } else {
         print  "NG1: $which2   $G1_L $G1_R   $G2_R $G2_L   $G1_B $G2_B   $G1_D $G2_D \n"; 
      return (1 == 0);        # false (doesn't look like a speciation)
   }
}




sub is_node_speciation{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $Sp_grp2 = shift;
   my $D1group = shift;
   my $D2group = shift;
   my $S12sum_range = shift;
   my $D12sum_range = shift;
print "S12sum range: $S12sum_range  D12sum range: $D12sum_range \n";
   #  my $group = shift;

print "s1s2d1d2:  $Sp_grp1  $Sp_grp2 $D1group $D2group \n";

   my $which = undef;
   my ($Lspecies_count, $Lcat_spcount, $cat_leafcountL, $Lid_presence,
       $Rspecies_count, $Rcat_spcount, $cat_leafcountR, $Rid_presence,
       $Bspecies_count, $Bcat_spcount, $cat_leafcountB, $Bid_presence)
     = $grouper->three_subtrees_species_counts($node, $sequenceid_species);
   my $Sp_grp1name = ($Sp_grp1 =~ /\s*([^,\s]+)/)? $1 : '';
   my $Sp_grp2name = ($Sp_grp2 =~ /\s*([^,\s]+)/)? $1 : '';
   my $D1groupname = ($D1group =~ /\s*([^,\s]+)/)? $1 : '';
   my $D2groupname = ($D2group =~ /\s*([^,\s]+)/)? $1 : '';
   my $species_count;
   #  my $outstring = '';
   my $seqid_presence;
   my %which_outstring = ();
   #   print "Top of is_node_speciation. ", scalar @{$node->get_implicit_names()}, "\n";

   #   my $Lgroupname = $Sp_grp1name;
   #   my $Rgroupname = $Sp_grp2name;
   #   my $Bgroupname = $D1groupname;

   my $L_S12sum_spcount = ($Lcat_spcount->{$Sp_grp1name} // 0)  +  ($Lcat_spcount->{$Sp_grp2name} // 0);
   my $R_S12sum_spcount = ($Rcat_spcount->{$Sp_grp1name} // 0)  +  ($Rcat_spcount->{$Sp_grp2name} // 0);
   my $B_S12sum_spcount = ($Bcat_spcount->{$Sp_grp1name} // 0)  +  ($Bcat_spcount->{$Sp_grp2name} // 0);

   my $L_D12sum_spcount = ($Lcat_spcount->{$D1groupname} // 0)  +  ($Lcat_spcount->{$D2groupname} // 0);
   my $R_D12sum_spcount = ($Rcat_spcount->{$D1groupname} // 0)  +  ($Rcat_spcount->{$D2groupname} // 0);
   my $B_D12sum_spcount = ($Bcat_spcount->{$D1groupname} // 0)  +  ($Bcat_spcount->{$D2groupname} // 0);

   my ($S12min_sp, $S12max_sp) = ($S12sum_range =~ /^\s*(\S+)\s*,\s*(\S+)\s*/)? ($1, $2) : (-1, 10000);
   my ($D12min_sp, $D12max_sp) = ($D12sum_range =~ /^\s*(\S+)\s*,\s*(\S+)\s*/)? ($1, $2) : (-1, 10000);

   my $LR_S12sum_spcount = $L_S12sum_spcount + $R_S12sum_spcount;
   my $RB_S12sum_spcount = $R_S12sum_spcount + $B_S12sum_spcount;
   my $BL_S12sum_spcount = $B_S12sum_spcount + $L_S12sum_spcount;

   my $LR_D12sum_spcount = $L_D12sum_spcount + $R_D12sum_spcount;
   my $RB_D12sum_spcount = $R_D12sum_spcount + $B_D12sum_spcount;
   my $BL_D12sum_spcount = $B_D12sum_spcount + $L_D12sum_spcount;

    print STDERR "S12 sp present; min, max, this subtree: $S12min_sp  $S12max_sp    $LR_S12sum_spcount     S1, S2 sp counts:  " ,
         ($Lcat_spcount->{$Sp_grp1name} // 0)  +  ($Rcat_spcount->{$Sp_grp1name} // 0), "   ",
        ($Lcat_spcount->{$Sp_grp2name} // 0)  +  ($Rcat_spcount->{$Sp_grp2name} // 0), "\n";
print STDERR "D12 min max this subtree:  $D12min_sp  $D12max_sp  $LR_D12sum_spcount \n";
   if (
       ($LR_S12sum_spcount >= $S12min_sp  and   $LR_S12sum_spcount <= $S12max_sp)
       and 
       ($LR_D12sum_spcount >= $D12min_sp  and   $LR_D12sum_spcount <= $D12max_sp)
       and
       (two_subtrees_speciation($Lcat_spcount, $Rcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "LR")) 
      ) {
      print STDERR  "LR  ", 
        $Lcat_spcount->{$Sp_grp1name} // 0, "  ", $Lcat_spcount->{$Sp_grp2name} // 0, "  ", $Lcat_spcount->{$D1groupname} // 0, "      ",
          $Rcat_spcount->{$Sp_grp1name} // 0, "  ", $Rcat_spcount->{$Sp_grp2name} // 0, "  ", $Rcat_spcount->{$D1groupname} // 0, "\n";
      #  print $node->recursive_subtree_newick(), "\n";
      $which = 'LR';
      $species_count = add_hashes($Lspecies_count, $Rspecies_count);
      $seqid_presence = add_hashes($Lid_presence, $Rid_presence);
      $grouper->group_species_counts($species_count);
      my $outstring = '';
      for my $gname ($Sp_grp1name, $Sp_grp2name) {
         $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
         my ($sumstr, $counts_str) = $grouper->group_speciescount_strings($gname);
         $outstring .= sprintf("%s     ", $counts_str); # join("  ", $grouper->group_speciescount_strings($gname)) );
         #        $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
      }
      #   $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);

      $outstring .= 
        # "  LR  " . 
        $node->recursive_subtree_newick() . "\n";
      $which_outstring{'LR'} = $outstring;
   }
   if ($LROK_descend_into_children or !defined $which) { # if $LROK_descend_into_children is false, and we found a LR speciation for this node, then don't do the following.

      if (
          ($RB_S12sum_spcount >= $S12min_sp  and   $RB_S12sum_spcount <= $S12max_sp)
          and 
          ($RB_D12sum_spcount >= $D12min_sp  and   $RB_D12sum_spcount <= $D12max_sp)
          and
          two_subtrees_speciation($Rcat_spcount, $Bcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "RB")
         ) {
         print STDERR  "RB  ",
           $Rcat_spcount->{$Sp_grp1name} // 0, "  ", $Rcat_spcount->{$Sp_grp2name} // 0, "  ", $Rcat_spcount->{$D1groupname} // 0, "      ",
             $Bcat_spcount->{$Sp_grp1name} // 0, "  ", $Bcat_spcount->{$Sp_grp2name} // 0, "  ", $Bcat_spcount->{$D1groupname} // 0, "\n";
         $which = 'RB';
         $species_count = add_hashes($Rspecies_count, $Bspecies_count);
         $seqid_presence = add_hashes($Rid_presence, $Bid_presence);
         $grouper->group_species_counts($species_count);
         my $outstring = '';
         for my $gname ($Sp_grp1name, $Sp_grp2name) {
            $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
            my ($sumstr, $counts_str) = $grouper->group_speciescount_strings($gname);
            $outstring .= sprintf("%s     ", $counts_str); # join("  ", $grouper->group_speciescount_strings($gname)) );
         }
         #    $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);

         my $the_tree_copy = $node->get_tree()->copy();
         my $implicit_name_string = join( "\t", @{ $node->get_implicit_names() } );
         my $node_in_copy = $the_tree_copy->node_from_implicit_name_string($implicit_name_string);
         #   print "node in copy name: ", $node_in_copy->get_name(), "\n";
         my @chs = $node_in_copy->get_children();
         my $new_outnode = $chs[0]; # ${$node->get_children()}[0];
         $the_tree_copy->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
         @chs = $the_tree_copy->get_root()->get_children();
         my $new_outgroup_root = $chs[0]; # $the_tree->get_root()->get_children()[0];
         my $new_subtree_root = $chs[1]; # $the_tree->get_root()->get_children()[1];

         $outstring .= 
           # " RB  " . 
           $new_subtree_root->recursive_subtree_newick() . "\n";
         $which_outstring{'RB'} = $outstring;
      }

      if (
          ($BL_S12sum_spcount >= $S12min_sp  and   $BL_S12sum_spcount <= $S12max_sp)
          and 
          ($BL_D12sum_spcount >= $D12min_sp  and   $BL_D12sum_spcount <= $D12max_sp)
          and
          two_subtrees_speciation($Bcat_spcount, $Lcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "BL")
         ) {
         print STDERR  "BL  ",
           $Bcat_spcount->{$Sp_grp1name} // 0, "  ", $Bcat_spcount->{$Sp_grp2name} // 0, "  ", $Bcat_spcount->{$D1groupname} // 0, "      ",
             $Lcat_spcount->{$Sp_grp1name} // 0, "  ", $Lcat_spcount->{$Sp_grp2name} // 0, "  ", $Lcat_spcount->{$D1groupname} // 0, "\n";
         $which = 'BL';
         $species_count = add_hashes($Bspecies_count, $Lspecies_count);
         $seqid_presence = add_hashes($Bid_presence, $Lid_presence);
         #     while(my($k,$v) = each %$species_count){ print "$k:$v  "; } print "\n";
         $grouper->group_species_counts($species_count);
         my $outstring = '';
         for my $gname ($Sp_grp1name, $Sp_grp2name) {
            $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
            my ($sumstr, $counts_str) = $grouper->group_speciescount_strings($gname);
            $outstring .= sprintf("%s     ", $counts_str); # join("  ", $grouper->group_speciescount_strings($gname)) );
            #        $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
         }
         #     $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);

         my $the_tree_copy = $node->get_tree()->copy();
         my $implicit_name_string = join( "\t", @{ $node->get_implicit_names() } );
         my $node_in_copy = $the_tree_copy->node_from_implicit_name_string($implicit_name_string);
         #   print "node in copy name: ", $node_in_copy->get_name(), "\n";
         my @chs = $node_in_copy->get_children();
         my $new_outnode = $chs[1]; # ${$node->get_children()}[0];
         $the_tree_copy->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
         @chs = $the_tree_copy->get_root()->get_children();
         my $new_outgroup_root = $chs[0]; # $the_tree->get_root()->get_children()[0];
         my $new_subtree_root = $chs[1]; # $the_tree->get_root()->get_children()[1];


         # my @chs = $node->get_children();
         # my $new_outnode = $chs[1]; # ${$node->get_children()}[0];
         # # my $new_outnode = $node->get_children()[0];
         # my $the_tree = $node->get_tree();
         # $the_tree->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
         # @chs =  $the_tree->get_root()->get_children();
         # my $new_outgroup_root = $chs[0]; # $the_tree->get_root()->get_children()[0];
         # my $new_subtree_root = $chs[1]; # $the_tree->get_root()->get_children()[1];
         $outstring .= 
           # " BL  " . 
           $new_subtree_root->recursive_subtree_newick() . "\n";

         # $outstring .= "  BL  " . $node->recursive_subtree_newick() . "\n";
         $which_outstring{'BL'} = $outstring;
      }
   }
   #  if(exists $which_outstring{'LR'}){
   #     print $node->recursive_subtree_newick(), "\n";
   #  }
   return (\%which_outstring);
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


sub is_node_all_one_group{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $D1group = shift;
   my ($species_count, $cat_spcount, $cat_leafcount, $id_presence) =
     $grouper->species_and_group_counts($node->get_implicit_names(), $sequenceid_species);
   #print "S1, D1: $Sp_grp1  $D1group  ", $cat_spcount->{$Sp_grp1} // 0, "  ", $cat_spcount->{$D1group} // 0, " \n";
   my $OK = (
             (($cat_spcount->{$D1group} // 0)  == 0)
             and 
             (($cat_spcount->{$Sp_grp1} // 0) > 0)
            );
   #print STDERR "In is_node_all_one_group. OK: ", $OK? 'OK' : 'NG', " \n";
   return $OK;
}

sub recursive_is_node_all_one_group{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $D1group = shift;
   my $tree_filename = shift;
   my $allowed_only_tree = shift;

   my $OK = is_node_all_one_group($node, $grouper, $sequenceid_species, $Sp_grp1, $D1group);
   if ($OK) {
      my $n_leaves = scalar @{$node->get_implicit_names()};
      print "$tree_filename   ", $allowed_only_tree? '1' : '0', "  $n_leaves    ", $node->recursive_subtree_newick(), "\n";
   } else {
      #    if (1 or (! ($is_spec and ($which eq 'LR')))) { # if child subtrees (i.e. LR) look like speciation don't descend to child nodes.
      my @children = $node->get_children();
      for (@children) {
         recursive_is_node_all_one_group($_, $grouper, $sequenceid_species, $Sp_grp1, $D1group, $tree_filename, $allowed_only_tree);
      }
   }
}


sub analyze_tree_4{
# reroot near disallowed leaf (if any present), then
# just look for maximal subtrees with none from D1group.
   # output number of species and sequences in S groups 
   my $tree = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift; # keys: sequence ids, values: species
   my $group_spcount_ranges = shift;  # e.g. [ [erysimum,14,1000], [brassicaceae7,0,5] ]
   my $group_leafcount_ranges = shift; # e.g. [ [erysimum,30,10000], [brassicaceae7,0,9] ]
   my $Dgroups = shift; # e.g. ['basals', non_ams]
   my $tree_filename = shift;
   print STDERR "xx Dgroups ", join(", ", @$Dgroups), "\n"; # exit;
#   print "Tree in original form:  ", $tree->get_root()->recursive_subtree_newick(), "\n";


# get list of disallowed leaves:
   my $disallowed_leaves = get_species_group_leaves( $tree, $grouper, $Dgroups); # $grouper->get_group($D1group) );
   my $n_disallowed_leaves = scalar @$disallowed_leaves;
   print STDERR "N disallowed leaves: $n_disallowed_leaves \n";
   my $allowed_only_tree = ($n_disallowed_leaves == 0);
   if ($n_disallowed_leaves > 0) { # reroot with a disallowed-species leaf as outgroup:
      my $new_outnode = $disallowed_leaves->[0];
      print STDERR "# rerooting with ", $new_outnode->get_name(), " as outgroup.\n";
      $tree->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
   } else {                     # midpoint reroot
      print STDERR "# rerooting midpoint.\n";
      $tree->reset_root_to_point_on_branch($tree->find_midpoint());
   }
   $tree->get_root()->recursive_implicit_names(); #
   $tree->get_root()->recursive_implicit_species(); #
#   print "possibly rerooted tree to analyze:  ", $tree->get_root()->recursive_subtree_newick(), "\n";

   my $node = $tree->get_root();
   recursive_is_node_all_allowed($node, $grouper, $sequenceid_species, 
                                 $group_spcount_ranges, $group_leafcount_ranges, $Dgroups, 
                                 $tree_filename, $allowed_only_tree); #, $Fgroupname);
   $tree->impose_branch_length_minimum(1);
   $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
}

sub is_node_all_allowed{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $group_spcount_ranges = shift;
   my $group_leafcount_ranges = shift;
   my $Dgroups = shift;
   my ($all_allowed, $too_few_sp, $too_many_sp, $too_few_leaves, $too_many_leaves) = 
     (1, 0, 0, 0, 0);
   my ($species_count, $cat_spcount, $cat_leafcount, $id_presence) =
     $grouper->species_and_group_counts($node->get_implicit_names(), $sequenceid_species);
   print STDERR "N species in subtree: ", scalar keys %$species_count, "\n";
   print STDERR "     ", join(" ", keys %$species_count), "\n";
   my $n_leaves = scalar  @{$node->get_implicit_names()};

   print STDERR "yy Dgroups ", join(", ", @$Dgroups), "\n"; # exit;
   for my $dgrp (@$Dgroups) {
      my $Dspcount = $cat_spcount->{$dgrp} // 0; # number of sp present of $dgrp
      print STDERR "asdf: $dgrp $Dspcount \n";
      if ($Dspcount > 0) {
         $all_allowed = 0;
       #  print STDERR "has disallowed:  $has_disallowed \n";
         last;
      }
   }
   if ($all_allowed) {
    #  ($too_few_sp, $too_many_sp, $too_few_leaves, $too_many_leaves) = (undef, undef, undef, undef);
#   } else {
   #   ($too_few_sp, $too_many_sp, $too_few_leaves, $too_many_leaves) = (0, 0, 0, 0);
      #    my $s_leafcount = 0;
      for my $g_lcr (@$group_leafcount_ranges) {
         my ($g, $min, $max) = @$g_lcr;
         print STDERR "grp, min max leaf counts: $g $min $max \n";
         my $g_leafcount = $cat_leafcount->{$g} // 0;
         #   $s_leafcount += $g_leafcount;
         $too_few_leaves = 1 if($g_leafcount < $min);
         $too_many_leaves = 1 if($g_leafcount > $max);
      }
      for my $grp_scr (@$group_spcount_ranges) {
         my ($g, $min, $max) = @$grp_scr;
         my $g_spcount = $cat_spcount->{$g} // 0;

         print STDERR "grp, min max sp counts: $g $min $max   $g_spcount \n";
         #   $s_leafcount += $g_leafcount;
         $too_few_sp = 1 if($g_spcount < $min);
         $too_many_sp = 1 if($g_spcount > $max);
      }
   }
      print STDERR "In is_node_all_allowed: $all_allowed, $too_few_sp, $too_many_sp, $too_few_leaves, $too_many_leaves \n";
      return ($all_allowed, $too_few_sp, $too_many_sp, $too_few_leaves, $too_many_leaves);
}

sub recursive_is_node_all_allowed{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $group_spcount_ranges = shift;
   my $group_leafcount_ranges = shift;
   my $Dgroups = shift;
   my $tree_filename = shift;
   my ($all_allowed, $too_few_sp, $too_many_sp, $too_few_leaves, $too_many_leaves) = is_node_all_allowed($node, $grouper, $sequenceid_species, 
                                                                       $group_spcount_ranges, $group_leafcount_ranges, $Dgroups);

   if ($too_few_sp or $too_few_leaves) { # subtree has nothing of interest. no output. do not descend further.
      return;
   } elsif ((! $all_allowed) or $too_many_sp or $too_many_leaves) { # no output, but do descend into child subtrees.
      my @children = $node->get_children();
      for (@children) {
         recursive_is_node_all_allowed($_, $grouper, $sequenceid_species, 
                                       $group_spcount_ranges, $group_leafcount_ranges, $Dgroups, $tree_filename);
      }
   } else {                     # (1, 0, 0, 0, 0) => OK
      my $n_leaves = scalar @{$node->get_implicit_names()};
      $grouper->group_species_counts($node->get_implicit_species());
      my $outstring = '';
      $outstring .= sprintf( "%s   %i   %i",$tree_filename, $n_leaves, scalar @$group_leafcount_ranges);
      my $s_species_present = 0;
      for my $ansgrp (@$group_leafcount_ranges) {
         my $sgrp = $ansgrp->[0];
         my $n_sp_in_group = $grouper->get_n_species_in_group($sgrp);
         my ($s_summary_str, $s_spcount_str) =  $grouper->group_speciescount_strings($sgrp);
         my @summary_numbers = split(" ", $s_summary_str);
         $s_species_present += $summary_numbers[0];
         $outstring .= sprintf("    %i   %s",$n_sp_in_group, $s_spcount_str); # if($summary_numbers[0] >= $min_s_species_present);
      }
      $outstring .= "   " . sprintf("%s", $node->recursive_subtree_newick());
      print "$outstring \n"; # if($s_species_present >= $min_s_species_present);
   }
}


sub A_both_B_counts{
   my $key_region = shift;
   my ($nA, $nB, $nAB, $nOther) = (0, 0, 0, 0);
   while (my($k, $v) = each %$key_region) {
      if ($v eq 'A') {
         $nA++;
      } elsif ($v eq 'B') {
         $nB++;
      } elsif ($v eq 'AB') {
         $nAB++;
      } else {
         $nOther++;
      }
   }
   return ($nA, $nAB, $nB, $nOther);
}



# **********************************************
# another analysis type to look for duplication followed by loss in
# some lineages.

sub is_node_dl{ # looks like duplication followed by loss of one group at this node,
   # i.e. D1group absent, Sp_grp1 pres (most species) in both child subtrees,
   # Sp_grp2 pr
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $Sp_grp2 = shift;
   my $D1group = shift;
   my $D2group = shift;
   my $Ogroupname = shift;
   #  my $group = shift;

   my $which = undef;
   my ($Lspecies_count, $Lcat_spcount, $Lid_presence, 
       $Rspecies_count, $Rcat_spcount, $Rid_presence,
       $Aspecies_count, $Acat_spcount, $Aid_presence) 
     = $grouper->three_subtrees_species_counts($node, $sequenceid_species);
   my $Lgroupname = ($Sp_grp1 =~ /\s*([^,\s]+)/)? $1 : '';
   my $Rgroupname = ($Sp_grp2 =~ /\s*([^,\s]+)/)? $1 : '';
   my $Bgroupname = ($D1group =~ /\s*([^,\s]+)/)? $1 : '';
   my $species_count;
   #  my $outstring = '';
   my $seqid_presence;
   my %which_outstring = ();
   #   print "Top of is_node_speciation. ", scalar @{$node->get_implicit_names()}, "\n";

   if (two_subtrees_speciation($Lcat_spcount, $Rcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "LR")) {
      print STDERR  "LR  ", 
        $Lcat_spcount->{$Lgroupname} // 0, "  ", $Lcat_spcount->{$Rgroupname} // 0, "  ", $Lcat_spcount->{$Bgroupname} // 0, "      ",
          $Rcat_spcount->{$Lgroupname} // 0, "  ", $Rcat_spcount->{$Rgroupname} // 0, "  ", $Rcat_spcount->{$Bgroupname} // 0, "\n";
      #  print $node->recursive_subtree_newick(), "\n";
      $which = 'LR';
      $species_count = add_hashes($Lspecies_count, $Rspecies_count);
      $seqid_presence = add_hashes($Lid_presence, $Rid_presence);
      $grouper->group_species_counts($species_count);
      my $outstring = '';
      for my $gname ($Lgroupname, $Rgroupname) {
         $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
         $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
      }
      $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);
      $which_outstring{'LR'} = $outstring;
   }
   if ($LROK_descend_into_children or !defined $which) { # if $LROK_descend_into_children is false, and we found a LR speciation for this node, then don't do the following.

      if (two_subtrees_speciation($Rcat_spcount, $Acat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "RB")) {
         print STDERR  "RB  ", 
           $Rcat_spcount->{$Lgroupname} // 0, "  ", $Rcat_spcount->{$Rgroupname} // 0, "  ", $Rcat_spcount->{$Bgroupname} // 0, "      ",
             $Acat_spcount->{$Lgroupname} // 0, "  ", $Acat_spcount->{$Rgroupname} // 0, "  ", $Acat_spcount->{$Bgroupname} // 0, "\n";
         $which = 'RB';
         $species_count = add_hashes($Rspecies_count, $Aspecies_count);
         $seqid_presence = add_hashes($Rid_presence, $Aid_presence);
         $grouper->group_species_counts($species_count);
         my $outstring = '';
         for my $gname ($Lgroupname, $Rgroupname) {
            $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
            $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
         }
         $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);
         $which_outstring{'RB'} = $outstring;
      }

      if (two_subtrees_speciation($Acat_spcount, $Lcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "BL")) {
         print STDERR  "BL  ",
           $Acat_spcount->{$Lgroupname} // 0, "  ", $Acat_spcount->{$Rgroupname} // 0, "  ", $Acat_spcount->{$Bgroupname} // 0, "      ",
             $Lcat_spcount->{$Lgroupname} // 0, "  ", $Lcat_spcount->{$Rgroupname} // 0, "  ", $Lcat_spcount->{$Bgroupname} // 0, "\n";
         $which = 'BL';
         $species_count = add_hashes($Aspecies_count, $Lspecies_count);
         $seqid_presence = add_hashes($Aid_presence, $Lid_presence);
         #     while(my($k,$v) = each %$species_count){ print "$k:$v  "; } print "\n";
         $grouper->group_species_counts($species_count);
         my $outstring = '';
         for my $gname ($Lgroupname, $Rgroupname) {
            $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
            $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
         }
         $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);
         $which_outstring{'BL'} = $outstring;
      }
   }
   #  if(exists $which_outstring{'LR'}){
   #     print $node->recursive_subtree_newick(), "\n";
   #  }
   return (\%which_outstring);
}




# *********************************

sub old_is_node_speciation{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $Sp_grp2 = shift;
   my $D1group = shift;
   my $D2group = shift;
   my $Ogroupname = shift;
   my $S12sum_range = shift;
   my $D12sum_range = shift;
   #  my $group = shift;

   my $which = undef;
   my ($Lspecies_count, $Lcat_spcount, $Lid_presence, 
       $Rspecies_count, $Rcat_spcount, $Rid_presence,
       $Bspecies_count, $Bcat_spcount, $Bid_presence) 
     = $grouper->three_subtrees_species_counts($node, $sequenceid_species);
   my $Lgroupname = ($Sp_grp1 =~ /\s*([^,\s]+)/)? $1 : '';
   my $Rgroupname = ($Sp_grp2 =~ /\s*([^,\s]+)/)? $1 : '';
   my $Bgroupname = ($D1group =~ /\s*([^,\s]+)/)? $1 : '';
   my $species_count;
   #  my $outstring = '';
   my $seqid_presence;
   my %which_outstring = ();
   #   print "Top of is_node_speciation. ", scalar @{$node->get_implicit_names()}, "\n";

   my $Sp_grp1name = $Lgroupname;
   my $Sp_grp2name = $Rgroupname;
   my $L_S12sum_spcount = ($Lcat_spcount->{$Sp_grp1name} // 0)  +  ($Lcat_spcount->{$Sp_grp2name} // 0);
   my $R_S12sum_spcount = ($Rcat_spcount->{$Sp_grp1name} // 0)  +  ($Rcat_spcount->{$Sp_grp2name} // 0);
   my $B_S12sum_spcount = ($Bcat_spcount->{$Sp_grp1name} // 0)  +  ($Bcat_spcount->{$Sp_grp2name} // 0);

   
   my ($S12min_sp, $S12max_sp) = ($S12sum_range =~ /^\s*(\S+)\s*,\s*(\S+)\s*/)? ($1, $2) : (-1, 10000);
   my ($D12min_sp, $D12max_sp) = ($D12sum_range =~ /^\s*(\S+)\s*,\s*(\S+)\s*/)? ($1, $2) : (-1, 10000);

   if ( (two_subtrees_speciation($Lcat_spcount, $Rcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "LR")) ) {
      print STDERR  "LR  ", 
        $Lcat_spcount->{$Lgroupname} // 0, "  ", $Lcat_spcount->{$Rgroupname} // 0, "  ", $Lcat_spcount->{$Bgroupname} // 0, "      ",
          $Rcat_spcount->{$Lgroupname} // 0, "  ", $Rcat_spcount->{$Rgroupname} // 0, "  ", $Rcat_spcount->{$Bgroupname} // 0, "\n";
      #  print $node->recursive_subtree_newick(), "\n";
      $which = 'LR';
      $species_count = add_hashes($Lspecies_count, $Rspecies_count);
      $seqid_presence = add_hashes($Lid_presence, $Rid_presence);
      $grouper->group_species_counts($species_count);
      my $outstring = '';
      for my $gname ($Lgroupname, $Rgroupname) {
         $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
         my ($sumstr, $counts_str) = $grouper->group_speciescount_strings($gname);
         $outstring .= sprintf("%s     ", $counts_str); # join("  ", $grouper->group_speciescount_strings($gname)) );
         #        $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
      }
      #   $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);

      $outstring .= 
        # "  LR  " . 
        $node->recursive_subtree_newick() . "\n";
      $which_outstring{'LR'} = $outstring;
   }
   if ($LROK_descend_into_children or !defined $which) { # if $LROK_descend_into_children is false, and we found a LR speciation for this node, then don't do the following.

      if (two_subtrees_speciation($Rcat_spcount, $Bcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "RB")) {
         print STDERR  "RB  ",
           $Rcat_spcount->{$Lgroupname} // 0, "  ", $Rcat_spcount->{$Rgroupname} // 0, "  ", $Rcat_spcount->{$Bgroupname} // 0, "      ",
             $Bcat_spcount->{$Lgroupname} // 0, "  ", $Bcat_spcount->{$Rgroupname} // 0, "  ", $Bcat_spcount->{$Bgroupname} // 0, "\n";
         $which = 'RB';
         $species_count = add_hashes($Rspecies_count, $Bspecies_count);
         $seqid_presence = add_hashes($Rid_presence, $Bid_presence);
         $grouper->group_species_counts($species_count);
         my $outstring = '';
         for my $gname ($Lgroupname, $Rgroupname) {
            $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
            my ($sumstr, $counts_str) = $grouper->group_speciescount_strings($gname);
            $outstring .= sprintf("%s     ", $counts_str); # join("  ", $grouper->group_speciescount_strings($gname)) );
         }
         #    $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);


         my $the_tree_copy = $node->get_tree()->copy();
         my $implicit_name_string = join( "\t", @{ $node->get_implicit_names() } );
         my $node_in_copy = $the_tree_copy->node_from_implicit_name_string($implicit_name_string);
         #   print "node in copy name: ", $node_in_copy->get_name(), "\n";
         my @chs = $node_in_copy->get_children();
         my $new_outnode = $chs[0]; # ${$node->get_children()}[0];
         $the_tree_copy->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
         @chs = $the_tree_copy->get_root()->get_children();
         my $new_outgroup_root = $chs[0]; # $the_tree->get_root()->get_children()[0];
         my $new_subtree_root = $chs[1]; # $the_tree->get_root()->get_children()[1];

         $outstring .= 
           # " RB  " . 
           $new_subtree_root->recursive_subtree_newick() . "\n";
         $which_outstring{'RB'} = $outstring;
      }

      if (two_subtrees_speciation($Bcat_spcount, $Lcat_spcount, $Sp_grp1, $Sp_grp2, $D1group, $D2group, "BL")) {
         print STDERR  "BL  ",
           $Bcat_spcount->{$Lgroupname} // 0, "  ", $Bcat_spcount->{$Rgroupname} // 0, "  ", $Bcat_spcount->{$Bgroupname} // 0, "      ",
             $Lcat_spcount->{$Lgroupname} // 0, "  ", $Lcat_spcount->{$Rgroupname} // 0, "  ", $Lcat_spcount->{$Bgroupname} // 0, "\n";
         $which = 'BL';
         $species_count = add_hashes($Bspecies_count, $Lspecies_count);
         $seqid_presence = add_hashes($Bid_presence, $Lid_presence);
         #     while(my($k,$v) = each %$species_count){ print "$k:$v  "; } print "\n";
         $grouper->group_species_counts($species_count);
         my $outstring = '';
         for my $gname ($Lgroupname, $Rgroupname) {
            $outstring .= sprintf("%1i  ", scalar $grouper->get_group($gname)->keys());
            my ($sumstr, $counts_str) = $grouper->group_speciescount_strings($gname);
            $outstring .= sprintf("%s     ", $counts_str); # join("  ", $grouper->group_speciescount_strings($gname)) );
            #        $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($gname)) );
         }
         #     $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $Ogroupname);

         my $the_tree_copy = $node->get_tree()->copy();
         my $implicit_name_string = join( "\t", @{ $node->get_implicit_names() } );
         my $node_in_copy = $the_tree_copy->node_from_implicit_name_string($implicit_name_string);
         #   print "node in copy name: ", $node_in_copy->get_name(), "\n";
         my @chs = $node_in_copy->get_children();
         my $new_outnode = $chs[1]; # ${$node->get_children()}[0];
         $the_tree_copy->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
         @chs = $the_tree_copy->get_root()->get_children();
         my $new_outgroup_root = $chs[0]; # $the_tree->get_root()->get_children()[0];
         my $new_subtree_root = $chs[1]; # $the_tree->get_root()->get_children()[1];


         # my @chs = $node->get_children();
         # my $new_outnode = $chs[1]; # ${$node->get_children()}[0];
         # # my $new_outnode = $node->get_children()[0];
         # my $the_tree = $node->get_tree();
         # $the_tree->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
         # @chs =  $the_tree->get_root()->get_children();
         # my $new_outgroup_root = $chs[0]; # $the_tree->get_root()->get_children()[0];
         # my $new_subtree_root = $chs[1]; # $the_tree->get_root()->get_children()[1];
         $outstring .= 
           # " BL  " . 
           $new_subtree_root->recursive_subtree_newick() . "\n";

         # $outstring .= "  BL  " . $node->recursive_subtree_newick() . "\n";
         $which_outstring{'BL'} = $outstring;
      }
   }
   #  if(exists $which_outstring{'LR'}){
   #     print $node->recursive_subtree_newick(), "\n";
   #  }
   return (\%which_outstring);
}


sub group_ranges{
   my $group_ranges_str = shift;
  my $group_ranges = [];

   print "group ranges string: [$group_ranges_str]\n";
  $group_ranges_str =~ s/\s+//g;    # remove whitespace
   my @group_count_ranges = split(";", $group_ranges_str);
   for my $group_countrange (@group_count_ranges) {
      my @cols = split(',', $group_countrange);
      my ($grp, $min, $max);
      if (scalar @cols >= 3) {
         ($grp, $min, $max) = @cols[0..2];
         $min = 0 if($min eq '-');
         $max = $bigger_than_any_tree if($max eq '-');
         warn "group_ranges_str: $group_ranges_str has unexpected format.\n" if(scalar @cols > 3);
      } elsif (scalar @cols == 2) {
         ($grp, $min) = @cols;
         $max = $bigger_than_any_tree;
      } elsif (scalar @cols == 1) {
         $grp = $cols[0];
         ($min, $max) = (0, $bigger_than_any_tree);
      } else {
         die "group_ranges_str: $group_ranges_str has unexpected format.\n";
      }
    #  $sgrps->{$grp} = [$min, $max];
      print STDERR "#  group min max:  $grp $min $max \n";
      push @$group_ranges, [$grp, $min, $max];
   }
   return $group_ranges;
}


# sub analyze_tree_5{ # just look for maximal subtrees with none from D1group.
#    # output number of species and sequences in S groups 
#    my $tree = shift;
#    my $grouper = shift;         # Grouper object
#    my $sequenceid_species = shift;
#    my $Sp_grp1 = shift;
#    my $Sp_grp2 = shift;
#    my $Dgroup = shift;
#    my $tree_filename = shift;

#    my $disallowed_leaves_to_analyze = get_species_group_leaves( $tree, $grouper->get_group($Dgroup) );
#    my $n_disallowed_leaves = scalar @$disallowed_leaves_to_analyze;
#    my $allowed_only_tree = ($n_disallowed_leaves == 0);
#    if ($n_disallowed_leaves > 0) { # reroot with a disallowed-species leaf as outgroup:
#       my $new_outnode = $disallowed_leaves_to_analyze->[0];
#       print STDERR "# rerooting with ", $new_outnode->get_name(), " as outgroup.\n";
#       $tree->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
#       $tree->get_root()->recursive_implicit_names();   #
#       $tree->get_root()->recursive_implicit_species(); #
#    }

#    my $node = $tree->get_root();
#    recursive_is_node_X($node, $grouper, $sequenceid_species, $Sp_grp1, $Sp_grp2, $Dgroup, $tree_filename, $allowed_only_tree); #, $Fgroupname);
#    $tree->impose_branch_length_minimum(1);
#    $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
# }

sub recursive_is_node_X{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;
   my $Sp_grp2 = shift;
   my $Dgroup = shift;
   my $tree_filename = shift;
   my $allowed_only_tree = shift;

   my ($OK, $descend_further, $info_string, $sp_vennregion) = is_node_X($node, $grouper, $sequenceid_species, $Sp_grp1, $Sp_grp2, $Dgroup);
   # print "$OK    $descend_further    $info_string \n";
   if ($OK) {
      #   my $n_leaves = scalar @{$node->get_implicit_names()};
      #   $grouper->group_species_counts($node->get_implicit_species());
      print "$tree_filename    ";
      print "$info_string   ";
      # for my $sgrp (@$Sgroups) {
      #    my $n_sp_in_group = $grouper->get_n_species_in_group($sgrp);
      #    my ($s_summary_str, $s_spcount_str) =  $grouper->group_speciescount_strings($sgrp);
      #    print  "$n_sp_in_group  $s_spcount_str   ";
      # }
      print  $node->recursive_subtree_newick(), "\n";
   } elsif ($descend_further) {
      #    if (1 or (! ($is_spec and ($which eq 'LR')))) { # if child subtrees (i.e. LR) look like speciation don't descend to child nodes.
      my @children = $node->get_children();
      for (@children) {
         recursive_is_node_X($_, $grouper, $sequenceid_species, $Sp_grp1, $Sp_grp2, $Dgroup, $tree_filename, $allowed_only_tree);
      }
   } else {
      # reject this subtree, and no need to descend into child subtrees.
   }
}

sub is_node_X{
   my $node = shift;
   my $grouper = shift;         # Grouper object
   my $sequenceid_species = shift;
   my $Sp_grp1 = shift;         # e.g. erysimum
   my $Sp_grp2 = shift;         # e.g. other_brassicaceae
   my $Dgroup = shift;          # basals

   my $S1min = 15;
   my $S2min = 1;
   my $childS1min = 13;
   my $S1_both_child_min = 11;
   my $S2_both_child_max = 1;
   my ($species_count, $cat_spcount, $cat_leafcount, $id_presence) =
     $grouper->species_and_group_counts($node->get_implicit_names(), $sequenceid_species);
   my $S1_spcount = $cat_spcount->{$Sp_grp1} // 0;
   my $S2_spcount = $cat_spcount->{$Sp_grp2} // 0;
   my $D_spcount = $cat_spcount->{$Dgroup} // 0;
   my $D_OK = ($D_spcount == 0);
   my $S_OK = ( ($S1_spcount >= $S1min) and ($S2_spcount >= $S2min));
   #   my $OK = ( $D_OK and $S_OK );
   my $descend_further = 1;
   #print "s1spcount, s2spcount, Dspcount:  $S1_spcount, $S2_spcount, $D_spcount \n";
   if ( $D_OK and $S_OK ) {
      #   print "ok so far ...  $S1_spcount, $S2_spcount, $D_spcount \n";
      # check L, R child subtrees
      my @children = $node->get_children();
      my $n_children = scalar @children;
      if ($n_children == 2) {
         my ($S1_Lspcount, $S1_Rspcount, $S1_and_count, $S1_xor_count, $S1_sp_vennregion) = $grouper->group_species_venn($Sp_grp1, $children[0]->get_implicit_species(), $children[1]->get_implicit_species());
         return (0, 1, "$S1_Lspcount, $S1_Rspcount, $Sp_grp1 $S1_and_count $S1_xor_count", $S1_sp_vennregion) if($S1_and_count < $S1_both_child_min);
         my ($S2_Lspcount, $S2_Rspcount, $S2_and_count, $S2_xor_count, $S2_sp_vennregion) = $grouper->group_species_venn($Sp_grp2, $children[0]->get_implicit_species(), $children[1]->get_implicit_species());
         return (0, 1, "$S2_Lspcount, $S2_Rspcount, $Sp_grp2 $S2_and_count $S2_xor_count", $S2_sp_vennregion) if($S2_and_count > $S2_both_child_max);
         return (1, 0, "$Sp_grp1 $S1_and_count $S1_xor_count   $Sp_grp2 $S2_and_count $S2_xor_count", $S1_sp_vennregion);
      } else {
         warn "in is_node_X. Unexpected number of children: $n_children \n"; 
      }
   } elsif ($S_OK) { # but not $D_OK. child subtrees need to be checked.
      return (0, 1, "subtree contains disallowed. Bad, and descend into child subtrees.\n", {});
   } else {                 # subtree (and its subtrees) are ruled out
      return (0, 0, "Subtree ruled out based on species counts: $S1_spcount $S2_spcount $D_spcount", {});
   }
   # my @ns = map($cat_leafcount->{$_} // 0, @$Sgroups);
   # return ($OK, $D1_spcount, @ns);
}

# sub analyze_tree_1{
#    my $tree = shift;
#    my $grouper = shift;           # Grouper object
#    my $select_group_name = shift; # if undef, use 0th group of grouper obj.
#    my $sequenceid_species = shift;

#    my $n_disallowed_limit = shift;
#    my $min_sgroup_seqs = shift;
#    my $min_max_sgroup_exp = shift;
#    my $sgroup_max_to_median = shift;
#    my $disallowed_limit = shift // 0;

#    my $done_ids = {};        # keys are ids which are already analyzed
#    my @output_strings = ();

#    $select_group_name //= $grouper->get_ith_groupname(0); # specifies group upon which selection will be based. Default is first group.
#    my $leaves_to_analyze = get_species_group_leaves( $tree, $grouper->get_group($select_group_name) );

#    my $count_leaves_done = 0;
#    my $count_leaves_analyzed = 0;
#    my $reroot_count = 0;
#    for my $the_leaf (@$leaves_to_analyze) { # loop over leaves of interest
#       $count_leaves_done++;
     
#       my $sequence_id = $the_leaf->get_name();
#       #   print STDERR "  $sequence_id \n";
#       next if($done_ids->{$sequence_id}); # the subtree containing this one has already been found - skip.
#       $reroot_count++;
#       #  print STDERR "     rerooting near $sequence_id. $reroot_count \n";
#       ######  Reroot near the leaf  ######
#       $count_leaves_analyzed++;
#       my $seqid_presence = {$sequence_id => 1}; # this hash holds the ids in the maximal part of tree containing query and only pgroup sequences.

#       $tree->reset_root_to_point_on_branch($the_leaf, 0.5*$the_leaf->get_branch_length());
#       $tree->get_root()->recursive_implicit_names(); #

#       # print $tree->get_root()->recursive_subtree_newick(), "\n\n";

#       my @children = $tree->get_root()->get_children();
#       die "Root node has ", scalar @children, " children; should have exactly 2.\n" if(scalar @children != 2);
#       my $next_node = $children[1];
#       my $the_subtree_string = "[" . $the_leaf->get_name() . " ";
#       my $subtree_string_other_part = '';
#       my $n_dis_so_far = 0;
     
#       while (1) {
#          @children = $next_node->get_children();
#          my $n_children = scalar @children;
#          #$Lspecies_count : keys: species, values: number of leaves in subtree of that species.
#          if ($n_children == 2) {
#             my ($Lchild, $Rchild) = @children;
#             my ($Lspecies_count, $Lcat_spcount, $Lcat_leafcount, $Lid_presence) =  $grouper->species_and_group_counts($Lchild->get_implicit_names(), $sequenceid_species);
#             my ($Rspecies_count, $Rcat_spcount, $Rcat_leafcount, $Rid_presence) =  $grouper->species_and_group_counts($Rchild->get_implicit_names(), $sequenceid_species);
#             my $L_n_dis = $Lcat_leafcount->{'disallowed'} // 0; # number of leaves of disallowed species found in L subtree
#             my $R_n_dis = $Rcat_leafcount->{'disallowed'} // 0; # number of leaves of disallowed species found in R subtree
#             my $L_n_leaves = scalar keys %$Lid_presence;
#             my $R_n_leaves = scalar keys %$Rid_presence;
#             # max subtree with <= $n_disallowed_limit disallowed leaves
#             if ($L_n_dis + $R_n_dis + $n_dis_so_far <= (($n_disallowed_limit == 1)? 2 : 0)) { # add both, and done
#                $seqid_presence = increment_hash($seqid_presence, $Lid_presence);
#                $seqid_presence = increment_hash($seqid_presence, $Rid_presence);
#                $subtree_string_other_part = $Lchild->get_implicit_names()->[0] . " ()";
#                last;
#             } elsif ($L_n_dis + $n_dis_so_far <= $n_disallowed_limit) { # add L, continue
#                $seqid_presence = increment_hash($seqid_presence, $Lid_presence);
#                $next_node = $Rchild;
#                $subtree_string_other_part = $Lchild->get_implicit_names()->[0];
#                $n_dis_so_far += $L_n_dis;
#             } elsif ($R_n_dis + $n_dis_so_far <= $n_disallowed_limit) { # add R, continue
#                $seqid_presence = increment_hash($seqid_presence, $Rid_presence);
#                $next_node = $Lchild;
#                $subtree_string_other_part = $Rchild->get_implicit_names()->[0];
#                $n_dis_so_far += $R_n_dis;
#             } else {            # add neither, done
#                $subtree_string_other_part .= " (" . $next_node->get_implicit_names()->[0]. ")";
#                #    my $subtree_copy = $next_node->copy_subtree();
#                #    print $next_node->recursive_subtree_newick();
#                last;
#             }
#          } elsif ($n_children == 0) { # to get here next_node must be a single nonAM leaf -- done
#             last;
#          } elsif ($n_children > 2) {
#             die "node: ", $next_node->get_name(), " has ", scalar @children, " children (not binary tree).\n"
#          } else {
#             die "node: ", $next_node->get_name(), " has ", scalar @children, " children. ???.\n"
#          }
#       }                         # end of while(1) loop

#       my $species_count = ();
#       for my $k (keys %$seqid_presence) {
#          $done_ids->{$k} = 1;
#          $species_count->{$sequenceid_species->{$k}}++;
#       }

#       $grouper->group_species_counts($species_count);
#       my ($nsppres, $med_nseq, $max_nseq, $total_nseq) = $grouper->get_group_summary_info( $select_group_name );
#       if ( ($total_nseq >= $min_sgroup_seqs) and
#            ($max_nseq >= $min_max_sgroup_exp) and
#            ($max_nseq >= $sgroup_max_to_median*$med_nseq)
#          ) {
#          my $outstring = '';
#          for my $a_grpname ($grouper->get_groupnames()) {
#             next if($a_grpname eq 'disallowed');
#             $outstring .= sprintf("%1i  ", scalar $grouper->get_group($a_grpname)->keys());
#             $outstring .= sprintf("%s     ", join("  ", $grouper->group_speciescount_strings($a_grpname)) );
#          }
#          $outstring .= "  " . $grouper->ids_present_string($sequenceid_species, $seqid_presence, $select_group_name);
#          push @output_strings, $outstring;
#       }
#    }                            # loop over leaves to analyze
#    print STDERR "     rerooted $reroot_count times.\n"; #exit;
#    $tree->impose_branch_length_minimum(1);
#    $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
#    return \@output_strings;
# }


# sub analyze_tree_3{ # just look for subtrees with as many from Sp_grp1 as possible, with none from D1group.
#    my $tree = shift;
#    my $grouper = shift;         # Grouper object
#    my $sequenceid_species = shift;
#    my $Sp_grp1 = shift;
#    my $D1group = shift;
#    my $tree_filename = shift; 

#    #print "D1 group: $D1group \n";
#    my $disallowed_leaves_to_analyze = get_species_group_leaves( $tree, $grouper->get_group($D1group) );
#    #my @leafnames = map($_->get_name(), @$leaves_to_analyze);
#    #  print "LTA: ", join("; ", @leafnames), "\n";
#    my $n_disallowed_leaves = scalar @$disallowed_leaves_to_analyze;
#    print "$n_disallowed_leaves \n";
#    my $allowed_only_tree = ($n_disallowed_leaves == 0);
#    if ($n_disallowed_leaves > 0) { # reroot with a disallowed-species leaf as outgroup:
#       my $new_outnode = $disallowed_leaves_to_analyze->[0];
#       $tree->reset_root_to_point_on_branch($new_outnode, 0.5*$new_outnode->get_branch_length());
#       $tree->get_root()->recursive_implicit_names(); #
#       #  print $tree->get_root()->recursive_subtree_newick(), "\n";
#    }

#    my $node = $tree->get_root();
#    recursive_is_node_all_one_group($node, $grouper, $sequenceid_species, $Sp_grp1, $D1group, $tree_filename, $allowed_only_tree); #, $Fgroupname);
#    $tree->impose_branch_length_minimum(1);
#    $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
# }
