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
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use TomfyMisc qw ( newick_genspid2idgensp  increment_hash  format_newick_species_info );
use Grouper;
use Hash::Ordered;
use Set::Scalar;

# rta.pl
#
# Usage example: rta.pl -in '*.newick' -grp1 AMs -grp2 nonAMs -other others

# Defaults:
my $input_pattern = undef;
my $gg_filename   = undef;
my $group1_species = '';
my $group2_species = '';
my $other_species = '';
my $min_max_grp1_sequences = 5; # only output if group1 has at least one species with >= this number of leaves.
my $expansion_relative_to_median = 4; # only output if group1 has at least one species with >= this number times more leaves than the (group1) median.
#my $disallowed_species = undef; # 'Arabidopsis_thaliana,Solanum_lycopersicum'; # '13_nonAM_angiosperms'; # '20_nons';

# Process long cl options
GetOptions(
	   'input_pattern=s' => \$input_pattern, # e.g. '*.newick'
           'grp1=s' => \$group1_species, # filename or string (e.g. 'Coffea_arabica,Coffea_canephora') specifying positive species.
           'grp2=s' => \$group2_species, #
           'other=s' => \$other_species,
           'min_max=i' => \$min_max_grp1_sequences,
           'exp_rel_to_median=f' => \$expansion_relative_to_median,
           'ggfile=s'    => \$gg_filename, # defines species-sequence association.
	  );


######  store the gene-genome information  ######
my ($seqid_speciesgg, $species_countgg) = (defined $gg_filename)? store_gg_info($gg_filename) : ({}, {});
my @sspecies = sort {$a cmp $b}  keys %$species_countgg; # sorted species


######   get the input filenames (should be newick format)  ######
die "input pattern is undefined.\n" if(!defined $input_pattern);
my $input_filenames_str = `ls $input_pattern`;
die "No files found matching pattern: $input_pattern \n" if($input_filenames_str eq '');
my @input_filenames = split(" ", $input_filenames_str);


######  get the groups of species  ######
my ($group1, $grp1_sequencecount) = get_species_HO($group1_species); # get_group_species($group1_species); # $Hash::Ordered->new();
my ($group2, $grp2_sequencecount) = get_species_HO($group2_species); # Hash::Ordered->new();
my ($other, $other_sequencecount) =  get_species_HO($other_species); # Hash::Ordered->new();

print "# group1: \n#  ", join(", ", $group1->keys()), "\n";
print "# group2: \n#  ", join(", ", $group2->keys()), "\n";
print "# other: \n#  ", join(", ", $other->keys()), "\n";
my $size_ratio = $grp1_sequencecount / $grp2_sequencecount;
print "# group1, group2, other sequence counts: $grp1_sequencecount  $grp2_sequencecount  $other_sequencecount .  grp1/grp2 size ratio: $size_ratio \n"; 

my $category_species = Hash::Ordered->new(
                                          'group1' => $group1,
                                          'group2' => $group2,
                                          'other' => $other,
                                         );
# check that group1, group2, other don't overlap:
my $intrsctn;
if ( ($intrsctn = intersection($group1, $group2)) ne 'empty' ) {
   die "Species $intrsctn present in both 'group1' and 'group2'. Goodbye\n";
} elsif (  ($intrsctn = intersection($group1, $other)) ne 'empty' ) {
   die "Species $intrsctn present in both 'group1' and 'other'. Goodbye\n";
} elsif (  ($intrsctn = intersection($group2,$other)) ne 'empty' ) {
   die "Species $intrsctn present in both 'group2' and 'other'. Goodbye\n";
}

# get species_category (category == group) hash
my $species_category = {};
for my $a_grp ( $category_species->keys() ) { # ('group1', 'group2', 'other') {
   my $a_set = $category_species->get($a_grp);
   for my $a_species ($a_set->keys()) {
      $species_category->{$a_species} = $a_grp;
      print "# $a_species  $a_grp \n";
   }
} print "# \n";

print "# Will analyze ", scalar @input_filenames, " files.\n";
for my $input_filename (@input_filenames) {
   my ($newick, $seqid_species);
   print STDERR "$input_filename\n";
   open my $fh_in, "<", "$input_filename" or die "Couldn't open $input_filename for reading.\n";
   my $newick_expression = '';
   while (<$fh_in>) {
      $newick_expression .= $_;
   }
   close $fh_in;
   $newick_expression =~ s/\s+//g; # remove whitespace
   $newick_expression =~ s/;$//;   # remove final ;

   my @leaf_ids = ();
   my $leaf_count = 0;
   if ($newick_expression) {
   ($newick, $seqid_species) =  format_newick_species_info($newick_expression, 1);
      my $the_tree = make_tree($newick);
      my $grouper = Grouper->new($category_species);

      # whole tree:
      $grouper->group_species_counts($the_tree->get_root());
      my ($nsppres, $med_nseq, $max_nseq, $total_nseq) = $grouper->get_group_summary_info('group1');
      if ( ($nsppres >= 1) and ($max_nseq >= $min_max_grp1_sequences)  and  ($max_nseq >= $expansion_relative_to_median * $med_nseq) ) {
         print "$input_filename    ";
         for ('group1', 'group2') {
            printf("%1i ", scalar $category_species->get($_)->keys() );
            print join("  ", $grouper->group_speciescount_strings($_)), "     ";
         }
         # print $grouper->get_group_sequences_string($the_tree->get_root(), 'group1', $seqid_species);
         print "\n";
      }

      $the_tree->impose_branch_length_minimum(1);
      $the_tree->decircularize(); # done with tree - decircularize so can be garbage collected.
      #  my $output_string = analyze_tree($the_tree, $group1, $group2, $seqid_species, $species_category);
   } else {                     # the newick string is empty
      print STDERR "File $input_filename does not contain a newick string.\n";
   }

}                               # end of loop over input files
print "# Done.\n";
# end of main


##############  Subroutines  ################

sub intersection{ # return 'empty' if intersection is empty, otherwise first value encountered which is in both.
   my $h1 = shift;              # Hash::Ordered
   my $h2 = shift;              # Hash::Ordered
   my @intersection_elements = ();
   for my $v ($h1->keys()) {
      push @intersection_elements, $v if(defined $h2->get($v));
   }
   if (scalar @intersection_elements  >  0) {
      return join(",", @intersection_elements);
   } else {
      return 'empty';
   }
}

sub get_species_group_leaves{ # return a arrayref of node objects for the leaves of species in group.
   my $tree = shift;
   my $group = shift;    # Hash::Ordered of species (values are all 1)
   my @selected_leaves = ();
   for my $a_leaf ($tree->get_leaves()) {
      my $leaf_species = $a_leaf->get_species();
      #  if (exists $group->{$leaf_species}) {
      if ($group->contains($leaf_species)) {
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

sub wholetree_species_content{
   my $tree = shift;

   my %species_count = ();
   my @species = @{$tree->get_root()->get_implicit_species()};
   for my $sp (@species) {
      $species_count{$sp}++;
   }
   return \%species_count;
}

sub analyze_tree{
   my $tree = shift;
   my $group1 = shift;  #
   my $group2 = shift;
   my $sequenceid_species = shift; # keys: sequence ids, values: species
   my $species_category = shift; # keys: species, values: categories (groups)
   #   my $sizeratio = shift; # number of sequences in group1/ number of sequences in group2
   my $done_ids = {};        # keys are ids which are already analyzed
   my $output_string = '';

   my $root = $tree->get_root();
   select_interesting_wholetree($tree, $species_category);

   my %species_count = ();
   my @species = @{$root->get_implicit_species()};
   for my $sp (@species) {
      $species_count{$sp}++;
   }
   # while(my($k,$v) = each %species_count){
   #    print "species, count:   $k $v \n";
   # }
   #  my ($wholetree_grp_spcount, $wholetree_grp_seqcount) = subtree_group_counts($root, $species_category);

   # $tree->postorder_traversal( \&select_interesting_subtree, $root, $species_category, $wholetree_grp_spcount, $wholetree_grp_seqcount );
}

sub subtree_group_counts{ # for the subtree below $node, get the number of species and sequences in each group.
   my $node = shift;
   my $species_cat = shift;
   #  my $size_ratio = shift;
   #   my $min_group1_seqcount = shift || 10;
   my %grp_seqcount = ('group1' => 0, 'group2' => 0, 'other' => 0);
   my %grp_spcount = ('group1' => 0, 'group2' => 0, 'other' => 0);
   my %species_count = ();
   for my $sp ( @{$node->get_implicit_species()} ) {
      $species_count{$sp}++;
      my $grp = $species_cat->{$sp} // 'unknown';
      if ($species_count{$sp} == 1) {
         $grp_spcount{$grp}++;
      }
      $grp_seqcount{$grp}++;
   }
   # select_interesting_subtree($node, 
   return (\%grp_spcount, \%grp_seqcount);
}

sub print_expansion_information{
   my $grp_spcount = shift;
   my $grp_seqcount = shift;
   my $grp1_expansion = shift; # $grp_seqcount->{group1}/$grp_spcount->{group1};
   my $grp2other_expansion = shift; # (($grp_spcount->{group2} + $grp_spcount->{other}) == 0)?
   my $subtree_str = shift;

   my $exp_ratio = ($grp2other_expansion > 0)?  $grp1_expansion/$grp2other_expansion : 1000;
   #    0 : ($grp_seqcount->{group2} + $grp_seqcount->{other}) / ($grp_spcount->{group2} + $grp_spcount->{other});
   for ('group1', 'group2', 'other', 'unknown') {
      printf("%2i ", $grp_spcount->{$_} // 0);
   }
   print "  ";
   for ('group1', 'group2', 'other', 'unknown') {
      printf("%2i ", $grp_seqcount->{$_} // 0);
   }
   print "  ";
   printf("%5.2f %5.2f %5.2f   %s\n", $grp1_expansion, $grp2other_expansion, $exp_ratio, $subtree_str);
}

  sub select_interesting_wholetree{
     my $tree = shift;
     my $species_cat = shift;
     my $min_group1_seqcount = shift // 8;
     my $min_expansion_ratio = shift // 3;

     my $whole_tree_leaf_count = $tree->get_leaf_count(); #$node->get_tree()->get_leaf_count();
     my ($grp_spcount, $grp_seqcount) = subtree_group_counts($tree->get_root(), $species_cat);

     if ($grp_seqcount->{group1} > $min_group1_seqcount) { # tree has enough from group1
        my ($nAsp, $nBsp) = ($grp_spcount->{group1},  $grp_spcount->{group2} + $grp_spcount->{other});
        my ($nAseq, $nBseq) = ($grp_seqcount->{group1},  $grp_seqcount->{group2} + $grp_seqcount->{other});
             print STDERR "$nAsp  $nBsp    $nAseq  $nBseq \n";
        if ( ($nBseq == 0)  or  ($nAseq/$nAsp >= $min_expansion_ratio * $nBseq/$nBsp) ) {
           my $grp2oth_expansion = ($nBseq == 0)? 0.1 : $nBseq/$nBsp;
           print_expansion_information($grp_spcount, $grp_seqcount, $nAseq/$nAsp, $grp2oth_expansion, 'whole tree');
        }
     }
  }

sub select_interesting_subtree{
   my $node = shift;
   my $species_cat = shift;
   my $wholetree_grp_spcount = shift;
   my $wholetree_grp_seqcount = shift;
   my $min_group1_seqcount = shift || 8;
   my $min_expansion_ratio = shift || 3;

   my $whole_tree_leaf_count = $node->get_tree()->get_leaf_count();
   my $subtree_leaf_count = scalar @{$node->get_implicit_names()};
   my ($grp_spcount, $grp_seqcount) = subtree_group_counts($node, $species_cat);
   #     print $grp_seqcount->{group1}, "  ",  $grp_seqcount->{group2}, "  ",  $grp_seqcount->{other}, "\n";
   #   return if($grp_seqcount->{group1} < $min_group1_seqcount);

   if ($subtree_leaf_count == $whole_tree_leaf_count  and  $grp_seqcount->{group1} > $min_group1_seqcount) { # whole tree
      my $grp1_expansion = $grp_seqcount->{group1}/$grp_spcount->{group1};
      my $grp2other_expansion = ($grp_seqcount->{group2}+$grp_seqcount->{other})/($grp_spcount->{group2} + $grp_spcount->{other}) // 0;

      if ( ($grp_seqcount->{group1} >= $min_group1_seqcount) and $grp1_expansion >= $min_expansion_ratio * $grp2other_expansion) {

         print_expansion_information($grp_spcount, $grp_seqcount, $grp1_expansion, $grp2other_expansion, 'whole tree');
      }
   } else {
      #  my ($grp1_expansion, $grp2other_expansion);
      my ($left_rep, $right_rep);
      # print "Node name:  ", $node->get_name(), "\n";
      if ($node->is_leaf()) {
         ($left_rep, $right_rep) =  ($node->get_name(), "---");
      } else {                  # not leaf
         my @children = $node->get_children();
         ($left_rep, $right_rep) = ($children[0]->get_implicit_names()->[0], $children[1]->get_implicit_names()->[0]);
      }
      #   print $grp_spcount->{group1}, "  ",  $grp_spcount->{group2}, "  ",  $grp_spcount->{other}, "\n";
      if ( 0 and $grp_spcount->{other} >= 1
           and ($grp_seqcount->{group1} >= $min_group1_seqcount)) {
         my $grp1_expansion = $grp_seqcount->{group1}/$grp_spcount->{group1};
         my $grp2other_expansion = (($grp_spcount->{group2} + $grp_spcount->{other}) == 0)?
           0 : ($grp_seqcount->{group2}+$grp_seqcount->{other})/($grp_spcount->{group2} + $grp_spcount->{other});
         if ($grp1_expansion >= 3 * $grp2other_expansion) {
            #      print "YYY:   $grp1_expansion  $grp2other_expansion \n";
            print_expansion_information($grp_spcount, $grp_seqcount, $grp1_expansion, $grp2other_expansion, "$left_rep  $right_rep");
         }
      }
      my $complement_grp1_spcount = $wholetree_grp_spcount->{group1} - $grp_spcount->{group1};
      my $complement_grp2other_spcount = $wholetree_grp_spcount->{group2} - $grp_spcount->{group2}  +  $wholetree_grp_spcount->{other} - $grp_spcount->{other};
      my $complement_grp1_seqcount = $wholetree_grp_seqcount->{group1} - $grp_seqcount->{group1};
      my $complement_grp2other_seqcount = $wholetree_grp_seqcount->{group2} - $grp_seqcount->{group2}  +  $wholetree_grp_seqcount->{other} - $grp_seqcount->{other};

      my $grp1_expansion = ($complement_grp1_spcount == 0)? 0 : $complement_grp1_seqcount/$complement_grp1_spcount;
      my $grp2other_expansion = ($complement_grp2other_spcount  == 0)?
        0 : ($complement_grp2other_seqcount / $complement_grp2other_spcount);

      if ( 0 and
           ($grp_seqcount->{group2} + $grp_seqcount->{other} > 0) and # don't do complement if subtree is all group1
           (! $node->get_parent()->is_root()) # don't do complement for children of root             and ($grp_spcount->{group2} + $grp_spcount->{other} > 0)
           and ($complement_grp1_seqcount > $min_group1_seqcount) 
           and  ($wholetree_grp_spcount->{other} - $grp_spcount->{other}) >= 1) {
         if ($grp1_expansion >= 3 * $grp2other_expansion) {
            #       print "ZZZ:  $grp1_expansion  $grp2other_expansion \n";
            print_expansion_information($grp_spcount, $grp_seqcount,  $grp1_expansion, $grp2other_expansion, "$left_rep  $right_rep   C");
         }
      }
   }
}

sub species_category_id_info_string{
   my $ids_in_subtree = shift ;
   my $seqid_sp = shift;
   my $group = shift;         # Hash::Ordered; keys are species names.
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
   $output_string .= sprintf("  %2i %2i  ", scalar keys %groupspecies_counts, scalar @group_ids); # number of species in group present in subtree
   my @sorted_psp_counts = map($groupspecies_counts{$_} // 0, $group->keys());
   for my $count (@sorted_psp_counts) {
      $output_string .= sprintf("%1i ", $count);
   }
   $output_string .= sprintf(" [%s]", join(";", sort @group_ids));
   return $output_string;
}

sub species_and_category_counts{
   # for a subtree, find:
   # number of leaves of each species ($species_leafcount)
   # number of species present of each category (
   # hash of ids present. (keys are ids in subtree, values are all 1)
   my $subtree_root = shift;
   my $id_species = shift;
   my $species_category = shift; # keys are species names (e.g. 'Aquilegia_coerulea'), values categories (e.g. '35_AM_angiosperms')
   my $species_leafcount = {}; # counts number of leaves of each species in subtree
   my $category_speciescount = {}; # counts number of species of each category which are present in subtree
   my %subtree_id_presence = ();   # keys: ids in subtree, values: 1
   my @ids = @{$subtree_root->get_implicit_names()};
   for my $an_id (@ids) {

      $subtree_id_presence{$an_id} = 1;

      my $species = 'unknown';
      if (exists $id_species->{$an_id}) {
         $species =  $id_species->{$an_id};
      } else {
         print STDERR "id: [$an_id]; species is unknown.\n"; exit;
      }

      my $category = $species_category->{$species} // 'other';
      $category_speciescount->{$category}++ if(! exists $species_leafcount->{$species});
      $species_leafcount->{$species}++;
   }
   return ($species_leafcount, $category_speciescount, \%subtree_id_presence);
}

sub taxonify_newick {
   my $newick        = shift;
   my %seqid_species = %{ shift @_ };
   $newick =~ s/\s+$//;         # remove final whitespace
   my $new_newick = $newick;
   $new_newick =~ s/ ([\(,]) \s* ([^,\):]+) \s* ([,\):])/$1$2\n$3/xg; # add newline after each leaf
   my @leaf_lines = split( "\n", $new_newick );
   my $last_bit = pop @leaf_lines;
   for (@leaf_lines) {

      if (/\[species=.*\]/) { # species info already there for this leaf - leave it alone
      } else {
         / [\(,] \s* ([^\(,\):]+) \s* $/x;
         my $seq_id = $1;
         if ( exists $seqid_species{$seq_id} ) {
            my $species = $seqid_species{$seq_id};
            $seq_id .= '[species=' . $species . ']';
         } else {
            warn "sequence $seq_id; no corresponding species found.\n";
         }
         s/ ([\(,]) \s* ([^\(,\):]+) \s* $/$1$seq_id/x;
      }
   }
   $new_newick = join( '', @leaf_lines ) . $last_bit;
   return $new_newick;
}

sub store_gg_info {	 #xx    # read in gene-genome association file
   # store in hash (keys: seqids, values: species). e.g. ('ATxxxx' => Arabidopsis_thaliana');
   my $gg_filename   = shift;
   my %seqid_species = ();
   my %species_count = ();
   if ( defined $gg_filename ) {
      if ( -f $gg_filename ) {
         open my $fh_gg, "<", "$gg_filename" or die "Couldn't open $gg_filename for reading.\n";
         while (<$fh_gg>) {
            my @cols = split( " ", $_ );
            my $species = shift @cols;
            $species =~ s/:$//; # remove final colon if present.
            for (@cols) {
               if ( exists $seqid_species{$_} ) {
                  warn "key $_ already stored with species: ",
                    $seqid_species{$_}, ", new species would be $species.\n";
               } else {
                  $seqid_species{$_} = $species;
                  $species_count{$species}++;
               }
            }
         }
         close $fh_gg; # done storing gg_file info in hash %seqid_species
      } else {         # 
         die "$gg_filename: no such file.\n";
      }
   } else {
      die "gg filename undefined. \n";
   }
   return (\%seqid_species, \%species_count);
}

# sub valuez{ # returns value of key/value pair, or 0 if pair doesn't exist for that key.
#    my $href = shift;
#    my $key = shift;
#    return (exists $href->{$key})? $href->{$key} : 0;
# }

sub get_group_species{
   my $group_species = shift; # filename, or e.g. 'genus_species1, genus_species2'
   my $the_group = Hash::Ordered->new(); # key: species; value: 1
   my $species_string = '';
   if (-f $group_species) {    # get the groups of species from a file
      open my $fh_grsp, "<", "$group_species";
      while (my $line = <$fh_grsp>) {
         $line =~ s/#.*$//;
         #  $line =~ s/,/ /g;
         #  for (split(" ", $line)) {
         #     $the_group->set($_ => 1);
         #  }
         $species_string .= $line;
      }
   } else {     # get the groups of species from command line argument
      $species_string = $group_species;
   }
   $species_string =~ s/,/ /g;
   for (split(" ", $species_string)) {
      $the_group->set($_ => 1);
   }
   return $the_group;
}

sub get_species_set{
   my $species_list = shift; # filename, or e.g. 'genus_species1, genus_species2', or 'genus_species1 36574, genus_species2 23365' (with genome sizes)
   my $sequence_count = 0;
   my $species_string = '';
   if (-f $species_list) {     # get the groups of species from a file
      open my $fhin, "<", "$species_list";
      while (my $line = <$fhin>) {
         $line =~ s/#.*$//;
         $line =~ s/,/ /;
         print STDERR "$line \n";
         if ($line =~ /^\s*(\S+)\s+(\d+)\s*$/ ) { # 
            $species_string .= "$1 ";
            $sequence_count += $2;
         } elsif ($line =~ /^\s*(\S+)\s*$/ ) {
            $species_string .= "$1 ";
            $sequence_count += 1;
         }
      }
      $species_string =~ s/,/ /g;
   } else { # $species_list is not a filename. get the groups of species from command line argument
      $species_string = $species_list;
      $species_string =~ s/,/ /g;
   }
   my $the_set = Set::Scalar->new(split(" ", $species_string));
   return ($the_set, $sequence_count);
}

sub get_species_HO{ # returns a Hash::Ordered
   my $species_str = shift; # filename, or e.g. 'genus_species1, genus_species2', or 'genus_species1 36574, genus_species2 23365' (with genome sizes)
   my $sequence_count = 0;
   my $species_string = '';
   my $species = Hash::Ordered->new();
   if (-f $species_str) {     # get the groups of species from a file
      open my $fhin, "<", "$species_str";
      while (my $line = <$fhin>) {
         $line =~ s/#.*$//;
         $line =~ s/,/ /;
         print STDERR "$line \n";
         if ($line =~ /^\s*(\S+)\s+(\d+)\s*$/ ) { # 
            $species->set($1 => $2);
            $sequence_count += $2;
         } elsif ($line =~ /^\s*(\S+)\s*$/ ) {
            $species->set($1 => 1);
            $sequence_count += 1;
         }
      }
   } else { # $species_list is not a filename. get the groups of species from command line argument
      $species_str =~ s/,/ /g;
      my @sps = split(" ", $species_str);
      for(@sps){
         $species->set($_ => 1);
         $sequence_count++;
      }
   }
   return ($species, $sequence_count);
}

sub get_subtree_species_count{
   my $node = shift;
   my @leaf_species = split(" ", $node->get_implicit_species());
   print STDERR join("; ", @leaf_species), "\n";
   my %species_count = ();
   for (@leaf_species) {
      $species_count{$_}++;
   }
   while (my($sp,$c) = each %species_count) {
      print STDERR "$sp $c\n";
   }
   return \%species_count;
}

sub node_species_multisets{
   my $node = shift;
   my $whole_tree_mset = shift // Multiset->new( $node->get_tree()->get_root()->get_implicit_species() );
   my $node_mset = Multiset->new( $node->get_implicit_species() );
   my $rest_of_tree_mset = $whole_tree_mset->minus( $node_mset );
   my @children = $node->get_children();
   my $L_mset = Multiset->new( $children[0]->get_implicit_species() );
   my $R_mset = Multiset->new( $children[1]->get_implicit_species() );

   return ($rest_of_tree_mset, $L_mset, $R_mset);
}
