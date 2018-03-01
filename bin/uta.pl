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
use Hash::Ordered;

# uta: unrooted tree analyzer
# analyze a set of trees in an unrooted fashion.
# the general idea is that we are interested in finding subtrees
# each of which contains at least one sequence belonging to a set of
# species, and no more than a specified number of sequences belonging to
# some other set of species.
# e.g. we could require at least one sequence from a coffee species
# (say C. arabica, C. canephora or C. eugenioides), and zero sequences of
# tomato or arabidosis.
# For each coffee sequence, root the tree near that sequence
# then work down through the tree, until reaching a node where both
# subtrees contain disallowed species. The branch above that node defines
# a bipartition of the tree, with one part the maximal subtree containing the
# starting leaf, and no disallowed species.
#
# specify a set of files each containing a single newick expression;
# specify a set of leaves to investigate (e.g. all coffee sequences)
# specify a set of disallowed species.

# Usage example:
# uta.pl  -gg ggfile  -in '*.newick'  -pgr pgroup  -dis dspecies

# Defaults:
my $input_pattern = undef;
my $gg_filename   = undef;
my $p_groups = 'pgroups'; # this is a default filename specifying the positive species.
my $disallowed_species = undef; # 'Arabidopsis_thaliana,Solanum_lycopersicum'; # '13_nonAM_angiosperms'; # '20_nons';

# Process long cl options
GetOptions(
	   'input_pattern=s' => \$input_pattern, # e.g. '*.newick'
           'pgroups=s' => \$p_groups, # filename or string (e.g. 'Coffea_arabica,Coffea_canephora') specifying positive species.
           'disallowed_species=s' => \$disallowed_species, #filename or string specifying positive species; default: all species not in $p_groups 
           'ggfile=s'    => \$gg_filename, # defines species-sequence association.
	  );

######   get the input filenames (should be newick format)  ######
die "input pattern is undefined.\n" if(!defined $input_pattern);
my $input_filenames_str = `ls $input_pattern`;
die "No files found matching pattern: $input_pattern \n" if($input_filenames_str eq '');
my @input_filenames = split(" ", $input_filenames_str);

######  store the gene-genome information  ######
my ($seqid_species, $species_count) = (defined $gg_filename)? store_gg_info($gg_filename) : (undef, {});
my @sspecies = sort {$a cmp $b}  keys %$species_count; # sorted species

######  get the groups of positive species  ######
my $group1 = Hash::Ordered->new();
my $group2 = Hash::Ordered->new();
my $other = Hash::Ordered->new();
my $the_group = $group1;
if (-f $p_groups) {            # get the groups of species from a file
   open my $fh_psp, "<", "$p_groups";
   while (my $line = <$fh_psp>) {
      $the_group = $group2 if($line =~ s/^\s*group2\s//i);
      $line =~ s/^\s*group1\s//i; #
      $line =~ s/#.*$//;
      $line =~ s/,/ /g;
      for (split(" ", $line)) {
         $the_group->set($_ => 1);
      }
   }
} else {        # get the groups of species from command line argument
   $the_group = $group2 if($p_groups =~ s/^\s*group2\s*//i); # if first nonwhitespace is 'group2' (case insensitive), delete and set $the_group to $group2
   $p_groups =~ s/,/ /g;
   for (split(" ", $p_groups)) {
      $the_group->set($_ => 1);
   }
}

######  get the set of disallowed species  ######
my $disallowed = Hash::Ordered->new();
if (defined $disallowed_species) {
   if (-f $disallowed_species) { # get disallowed species from file
      open my $fh_nsp, "<", "$disallowed_species";
      while (my $line = <$fh_nsp>) {
         $line =~ s/#.*$//;
         $line =~ s/,/ /g;
         for (split(" ", $line)) {
            $disallowed->set($_ => 1);
         }
      }
   } else {                 # get disallowed species from command line
      $disallowed_species =~ s/,/ /g;
      for (split(" ", $disallowed_species)) {
         $disallowed->set($_ => 1);
      }
   }
} else {        # use default for disallowed - complement of %pspecies
   for my $sp (@sspecies) {
      if ((!defined $group1->get($sp))  and  (!defined $group2->get($sp))) {
         $disallowed->set($sp => 1);
      }
   }
}

my $category_species = Hash::Ordered->new(
                                          'group1' => $group1,
                                          'group2' => $group2,
                                          'other' => $other,
                                          'disallowed' => $disallowed
                                         );
# check that psp, nsp don't overlap (if they do, warn and remove from psp)
# check that all psp, nsp are in $species_count ( from gg ).

my $species_category = {};
for my $a_cat ($category_species->keys()) {
   for my $a_species ($category_species->get($a_cat)->keys()) {
      $species_category->{$a_species} = $a_cat;
      print "# $a_species  $a_cat \n";
   }
}

for my $a_cat ($category_species->keys()) {
   my $species = $category_species->get($a_cat);
}

print STDERR "# Will analyze ", scalar @input_filenames, " files.\n";
for my $input_filename (@input_filenames) {
   print STDERR "# Starting analysis of $input_filename.\n";
   open my $fh_in, "<", "$input_filename" or die "Couldn't open $input_filename for reading.\n";
   my $the_input_newick = '';
   while (<$fh_in>) {
      $the_input_newick .= $_;
   }
   $the_input_newick =~ s/\s//g;
   $the_input_newick =~ s/;?\s*$//; # remove final ; and whitespace if present

   if ($the_input_newick) {
      if ($the_input_newick =~ /\[species=[a-zA-Z_]+\]/) {
         my $wrk = $the_input_newick;
         while ($wrk =~ s/([(,])([^(,:[]+)\[species=([^\]]+)\][:]/$1:/) {
            $seqid_species->{$2} = $3;
           # print "AAA: $wrk \n";
         #   print "id: $2,  species: $3 \n";
         }
      } else {
         if (scalar keys %$seqid_species > 0) { # use info from ggfile
            $the_input_newick = taxonify_newick( $the_input_newick, $seqid_species );
         } else {
         #   print $the_input_newick, "\n";
            $the_input_newick =~ s/([(,])(\S+)__([^:]+)/$1 $2 [species=$3]/g; # from id __genus_species: to id[species=genus_species]:
            $seqid_species->{$2} = $3;
            $the_input_newick =~ s/\s+//g;
         #   print $the_input_newick, "\n";
         }
         
      }
   #   print "XXXXX: $the_input_newick \n";
      my $the_tree = make_tree($the_input_newick);
      my $output_string = analyze_tree($the_tree, $group1, $group2, $seqid_species, $species_category);
      print "$input_filename\n", "$output_string\n";
   } else {                     # the newick string is empty
      print STDERR "File $input_filename does not contain a newick string.\n";
   }
}                            # end of loop over input files
print STDERR "# Done.\n";
# end of main


##############  Subroutines  ################

sub get_species_group_leaves{
   my $tree = shift;
   my $group = shift;    # Hash::Ordered of species (values are all 1)
   my @selected_leaves = ();
   for my $a_leaf ($tree->get_leaves()) {
      my $leaf_species = $a_leaf->get_species();
      #  if (exists $group->{$leaf_species}) {
      if (defined $group->get($leaf_species)) {
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
   my $group1 = shift;          # Hash::Ordered of a group of species
   my $group2 = shift;
   my $sequenceid_species = shift;
   my $species_category = shift;
   my $done_ids = {};        # keys are ids which are already analyzed
   my $output_string = '';

   my $leaves_to_analyze = get_species_group_leaves($tree, $group1);
   push @$leaves_to_analyze, @{get_species_group_leaves($tree, $group2)};
   for my $the_leaf (@$leaves_to_analyze) { # loop over leaves of interest

      ######  Reroot near the leaf  ######
      my $sequence_id = $the_leaf->get_name();
      next if($done_ids->{$sequence_id}); # the subtree containing this one has already been found - skip.
      my $seqid_presence = {$sequence_id => 1}; # this hash holds the ids in the maximal part of tree containing query and only AM hosts.

      $tree->reset_root_to_point_on_branch($the_leaf, 0.5*$the_leaf->get_branch_length());
      $tree->get_root()->recursive_implicit_names(); # 

      my @children = $tree->get_root()->get_children();

      my $next_node = $children[1];
      while (1) {
         @children = $next_node->get_children();
         my $n_children = scalar @children;

         if ($n_children == 2) {
            my ($Lspecies_count, $Lcat_spcount, $Lid_presence) =  species_and_category_counts($children[0], $sequenceid_species, $species_category);
            my ($Rspecies_count, $Rcat_spcount, $Rid_presence) =  species_and_category_counts($children[1], $sequenceid_species, $species_category);
            my $L_OK = valuez($Lcat_spcount, 'disallowed') == 0; # number of disallowed species found in L subtree is zero
            my $R_OK = valuez($Rcat_spcount, 'disallowed') == 0; # number of disallowed species found in R subtree is zero
            if ($L_OK) {
               $seqid_presence = add_hashes($seqid_presence, $Lid_presence);
               if ($R_OK) { # both subtrees have no negatives (whole tree has no negatives)
                  $seqid_presence = add_hashes($seqid_presence, $Rid_presence);
                  last;
               } else {         # R subtree has negatives
                  $next_node = $children[1];
               }
            } else {            # L subtree has negatives
               if ($R_OK) {
                  $seqid_presence = add_hashes($seqid_presence, $Rid_presence);
                  $next_node = $children[0];
               } else {         # both subtrees have negatives -- done
                  last;
               }
            }
         } elsif ($n_children == 0) { # to get here next_node must be a single nonAM leaf -- done
            last;
         } elsif ($n_children > 2) {
            die "node: ", $next_node->get_name(), " has ", scalar @children, " children (not binary tree).\n"
         } else {
            die "node: ", $next_node->get_name(), " has ", scalar @children, " children. ???.\n"
         }
      }                         # end of while(1) loop

      for my $k (keys %$seqid_presence) {
         $done_ids->{$k} = 1;
      }

      my @ids_in_subtree = keys %$seqid_presence;
      $output_string .= species_category_id_info_string(\@ids_in_subtree, $sequenceid_species, $group1);
      if (scalar $group2->keys() > 0) {
         $output_string .= species_category_id_info_string(\@ids_in_subtree, $sequenceid_species, $group2);
      }
      $output_string .= "\n";
   }                            # loop over leaves to analyze
   $tree->impose_branch_length_minimum(1);
   $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
   return $output_string;
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
            $species =~ s/:$//;	# remove final colon if present.
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

sub valuez{ # returns value of key/value pair, or 0 if pair doesn't exist for that key.
   my $href = shift;
   my $key = shift;
   return (exists $href->{$key})? $href->{$key} : 0;
}

sub add_hashes{
   my $sc1 = shift;
   my $sc2 = shift;
   while ( my ($sp, $c) = each %$sc2) {
      $sc1->{$sp} += $c;
   }
   return $sc1;
}

