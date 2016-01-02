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

BEGIN {
   $bindir = dirname( abs_path(__FILE__) ); # this has to go in Begin block so happens at compile time
   $libdir = $bindir . '/../lib';
}
use lib $libdir;


use Getopt::Long;

use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;
use CXGN::Phylo::CladeSpecifier;

my $default_gg_file_path           = undef;
my $default_species_tree_file_path = undef;

# Usage example:
# clades.pl -input x.newicks              (will use defaults)
# clades.pl -input x.newicks -clade_spec '4monocots,3'   (using default for disallowed species)
# clades.pl -input x.newicks -clade_spec (Brachypodium_distachyon,Oryza_sativa),1  -disallowed 5brassicas
# clades.pl -input x.newicks -clade_spec (Brachypodium_distachyon,Oryza_sativa),1  -disallowed Brassica_rapa,Arabidopsis_thaliana
#
# Defaults:
my $input_newicks_filename;
my $gg_filename   = $default_gg_file_path;
my $reroot_method = undef; # default is do not reroot ( mindl rerooting is default for nj_ml_bs.pl )

my $clade_specifiers = 
  '23_AMp_dicots,10 : 22_AMp_dicots,10 : 9_monocots,4 : 7_basals,1 : 11_AMnegatives,1';
#'12pdicots,6 : 6monocots,3 : 8basals,1 : 9negatives,1'; 
# '8dicots,7 : 7dicots,6 : 4monocots,3 : Selaginella_moellendorffii,1'; # default clade specifiers.

# $clade_specifiers = '8dicots,7 : 7dicots,6 : 4monocots,1 : 4monocots,2 : 4monocots,3 : Selaginella_moellendorffii,1';
# This would mean we will look for a clade with 7 out of 8 dicots, for one with 6 out of 7 dicots, 
# for one with 3 monocots, and for one with Selaginella

my $category_species =
  { # hashref. keys are names of predef taxon groups; values are hashrefs (keys taxa, values 1)

   # Amborella:
   'Amborella' => {
                   'Amborella_trichopoda' => 1, 
                  },
   '12_AM_monocots' => {
                        # 12 AM monocots
                        'Elaeis_guineensis' => 1,
                        'Phoenix_dactylifera'     => 1, # date palm
                        'Musa_acuminata' => 1,          # banana
                        'Zea_mays'                => 1, # maize
                        'Panicum_virgatum' => 1,        # switchgrass
                        'Phyllostachys_heterocycla' => 1, # bamboo, AM ??
                        'Brachypodium_distachyon' => 1,
                        'Sorghum_bicolor'         => 1,
                        'Oryza_sativa'            => 1, # rice
                        'Hordeum_vulgare'         => 1, # barley
                        'Triticum_urartu'       => 1,   # 
                        'Panicum_hallii' => 1,
                       },
   '23_AM_dicots' => {
                      # 22 AM dicots
                      'Aquilegia_coerulea' => 1, # columbine

                      'Solanum_lycopersicum' => 1, # tomato
                      'Solanum_tuberosum'    => 1, # potato
                      'Mimulus_guttatus' => 1,     # monkeyflower
                      'Fraxinus_excelsior' => 1,   # Ash
                      'Sesamum_indicum' => 1,

                      'Vitis_vinifera'       => 1, # grape

                      'Glycine_max'          => 1, # soy
                      'Phaseolus_vulgaris' => 1,
                      'Lotus_japonicus' => 1,
                      'Medicago_truncatula' => 1,

                      'Populus_trichocarpa'  => 1, # poplar
                      'Ricinus_communis'     => 1, # castor
                      'Cucumis_sativus'      => 1, # cucumber
                      'Manihot_esculenta' => 1,
                      'Salix_purpurea' => 1,

                      'Theobroma_cacao' => 1,
                      'Carica_papaya' => 1,
                      'Eucalyptus_grandis' => 1,
                      'Gossypium_raimondii' => 1,
                      'Citrus_clementina' => 1,
                      'Citrus_sinensis' => 1,
                      # Lupine:
                      'Lupinus_angustifolius' => 1,
                     },
   # '35_AM_angiosperms' => {

   #                         # Amborella:
   #                         'Amborella_trichopoda' => 1, 
   #                         # 12 AM monocots
   #                         'Elaeis_guineensis' => 1,
   #                         'Phoenix_dactylifera'     => 1, # date palm
   #                         'Musa_acuminata' => 1,	   # banana
   #                         'Zea_mays'                => 1, # maize
   #                         'Panicum_virgatum' => 1, # switchgrass
   #                         'Phyllostachys_heterocycla' => 1, # bamboo, AM ??
   #                         'Brachypodium_distachyon' => 1,
   #                         'Sorghum_bicolor'         => 1,
   #                         'Oryza_sativa'            => 1, # rice
   #                         'Hordeum_vulgare'         => 1, # barley
   #                         'Triticum_urartu'       => 1,   # 
   #                         'Panicum_hallii' => 1,

   #                         # 22 AM dicots
   #                         'Aquilegia_coerulea' => 1, # columbine

   #                         'Solanum_lycopersicum' => 1, # tomato
   #                         'Solanum_tuberosum'    => 1, # potato
   #                         'Mimulus_guttatus' => 1,     # monkeyflower
   #                         'Fraxinus_excelsior' => 1,   # Ash
   #                         'Sesamum_indicum' => 1,

   #                         'Vitis_vinifera'       => 1, # grape

   #                         'Glycine_max'          => 1, # soy
   #                         'Phaseolus_vulgaris' => 1,
   #                         'Lotus_japonicus' => 1,
   #                         'Medicago_truncatula' => 1,

   #                         'Populus_trichocarpa'  => 1, # poplar
   #                         'Ricinus_communis'     => 1, # castor
   #                         'Cucumis_sativus'      => 1, # cucumber
   #                         'Manihot_esculenta' => 1,
   #                         'Salix_purpurea' => 1,

   #                         'Theobroma_cacao' => 1,
   #                         'Carica_papaya' => 1,
   #                         'Eucalyptus_grandis' => 1,
   #                         'Gossypium_raimondii' => 1,
   #                         'Citrus_clementina' => 1,
   #                         'Citrus_sinensis' => 1,
   #                         # Lupine:
   #                         'Lupinus_angustifolius' => 1,
   #                        },
   '7_nonangiosperms' => { # which branch off before the monocot-dicot split;  non-angiosperms.
                          'Ostreococcus_tauri' => 1,
                          'Ostreococcus_lucimarinus' => 1,
                          'Volvox_carteri' => 1,
                          'Chlamydomonas_reinhardtii'  => 1,
                          'Physcomitrella_patens'      => 1, 
                          'Selaginella_moellendorffii' => 1, 
                          'Picea_abies' => 1, # norway spruce
                         },
   '13_nonAM_angiosperms' =>  {
                               Brassica_rapa           => 1, # turnip
                               Arabidopsis_thaliana    => 1,
                               Arabidopsis_lyrata      => 1,
                               Thellungiella_halophila => 1,
                               Eutrema_salsugineum => 1,
                               Capsella_rubella        => 1,
                               'Tarenaya_hassleriana' => 1,
                               Beta_vulgaris => 1,
                               Spinacia_oleracea => 1,
                               Nelumbo_nucifera => 1,
                               Utricularia_gibba => 1,
                               'Dianthus_caryophyllus' => 1,
                               'Spirodela_polyrhiza' => 1, # duckweed - monocot
                               'Zostera_marina' => 1,
                              },
   # '20_nons' => {
   #               'Ostreococcus_tauri' => 1,
   #               'Ostreococcus_lucimarinus' => 1,
   #               'Volvox_carteri' => 1,
   #               'Chlamydomonas_reinhardtii'  => 1,
   #               'Physcomitrella_patens'      => 1, 
   #               'Selaginella_moellendorffii' => 1, 
   #               'Picea_abies' => 1,             # norway spruce
   #               'Brassica_rapa'           => 1, # turnip
   #               'Arabidopsis_thaliana'    => 1,
   #               'Arabidopsis_lyrata'      => 1,
   #               'Thellungiella_halophila' => 1,
   #               'Eutrema_salsugineum' => 1,
   #               'Capsella_rubella'        => 1,
   #               'Tarenaya_hassleriana' => 1,
   #               'Beta_vulgaris' => 1,
   #               'Spinacia_oleracea' => 1,
   #               'Nelumbo_nucifera' => 1,
   #               'Utricularia_gibba' => 1,
   #               'Dianthus_caryophyllus' => 1,
   #               'Spirodela_polyrhiza' => 1, # duckweed - monocot
   #               'Zostera_marina' => 1,
   #              }
  };

my $species_category = {};
for my $a_cat (keys %$category_species) {
   for my $a_spec (keys %{$category_species->{$a_cat}}) {
      $species_category->{$a_spec} = $a_cat;
      #    print "$a_spec  $a_cat \n";
   }
}
#exit;

# For each clade, count the number of these species present: (
my $disallowed_species       = '11_AMnegatives'; # '5brassicas';
my $max_nbs                  = 1000000;	# big number - do all bs replicates by default.
my $species_tree_newick_file = $default_species_tree_file_path;	# if undefined or
my $id_to_show = undef;
my $type_to_show = 'ALL';	# other possibilities: 'NJ', 'ML'
my $species_to_show = undef; # for each clade, output the ids of this species in the clade.
my $show_clade_species = 0;
my $max_allowed_nonAM = 0;
my $max_allowed_nonang = 1000;
my $disallowed_group = '20_nons';
# Process long cl options
GetOptions(
	   'input_newicks=s' => \$input_newicks_filename,
	   'clade_specifiers=s' => \$clade_specifiers,
           'disallowed_group=s' => \$disallowed_group, # def: '20_nons', or, e.g.: '13_nonAM_angiosperms' ...
	   'disallowed_species=s' =>
	   \$disallowed_species, # either: name of predefined taxon group, or comma separated taxa names:
	   # e.g. 'Arabidopsis_thaliana, Brassica_rapa'
	   'maxnbs=i' => \$max_nbs,
	   # The following 3 options typically will not be used because the gene-genome association and tree rooting has already
	   # been done at an earlier stage. One can however choose to reroot the trees in a different way at this point.
	   'ggfile=s'        => \$gg_filename, # defines species-sequence association.
	   # gg file has 1 line per species: species name followed by whitespace-separated sequence ids.
	   'reroot=s' => \$reroot_method, # selects rerooting method. options: none, midpoint, minvar, mindl. 
	   'speciestreefile=s'  => \$species_tree_newick_file, # to override built-in species tree with 52 species
	   # (see sub reroot below to see this default species tree). Species tree only needed if rerooting, and using
	   # the mindl method.
	   'id_to_show=s' => \$id_to_show, 
	   'type=s' => \$type_to_show,
	   'species_to_show=s' => \$species_to_show,
	   'show_clade_species!' => \$show_clade_species,
	  );

if ( defined $species_tree_newick_file and !-f $species_tree_newick_file ) {
   die "$species_tree_newick_file is not a regular file. Exiting\n";
   $species_tree_newick_file = undef;
}
my $species_tree_obj   = get_species_tree($species_tree_newick_file);

###################### print clade specifiers, disallowed species. ######################
# print "# $clade_specifiers \n#\n";
# $clade_specifiers =~ s/\s+//g;
# my @clade_specs = split( ":", $clade_specifiers );

# my @clade_spec_objs = ();
# for (@clade_specs) {
#    push @clade_spec_objs, CXGN::Phylo::CladeSpecifier->new( $_, $predefined_taxon_groups );
# }
# # while ( my ( $i, $cso ) = each @clade_spec_objs ) {
# for (my $i=0; $i < scalar @clade_spec_objs; $i++ ) {
#    my $cso = $clade_spec_objs[$i];
#    print "# Clade ", $i + 1, " specs: \n", "# ", $cso->as_string(), "#\n";
# }
# print "# Disallowed species: ", join( ", ", @disallowed_species ), "\n#\n";

# print "# For each clade 5 numbers are given: \n",
#   "#   query-claderoot distance; \n", 
#   "# query-claderoot 'epsilon product' (product of 1-branch_support); \n",
#   "#   N query-species seqs in clade; \n",
#   "#   N disallowed species in clade; \n",
#   "#   N leaves in clade. \n";

##########################################################################################

if ( !defined $input_newicks_filename or !-f $input_newicks_filename ) {
   die "No input newicks file specified, or file [$input_newicks_filename] does not exist.";
}
open my $fh_in, "<", "$input_newicks_filename"
  || die "Couldnt open $input_newicks_filename for reading.\n";
my ($seqid_species, $species_count) = store_gg_info($gg_filename);
my @sspecies = sort {$a cmp $b}  keys %$species_count; # sorted species
print "# ", scalar @sspecies, " species: \n# ", join("\n# ", @sspecies), "\n";
#exit;
my $BS_count = 0;
while (<$fh_in>) {
   if (/^Id\s+(\S+)\s.*fam_size:\s+(\d+)/) {
      my $query_sequence_id = $1; # query id
      my $family_size = $2;
      next if(defined $id_to_show and ($query_sequence_id ne $id_to_show));	
      my $next_line   = <$fh_in>;
      my $type        = undef;

      #print STDERR "NEXT LINE: [", $next_line, "]\n";
      if ( $next_line =~ /^((\S+)\s+)?\s*\(/ ) { # newick string on this line
         my $the_input_newick;
         if ( $next_line =~ /^\s*\(/ ) { # newick string on this line, without 'NJ', 'ML' etc. (it is ML)
            $the_input_newick = $next_line;
            $type             = 'ML';
         } elsif ( $next_line =~ /^\s*(\S+)\s+(\(.*)$/ ) {
            $type = $1;
            next if($type_to_show ne 'ALL' and $type ne $type_to_show);
            if ( $type eq 'NJ_BS' ) {
               $BS_count++;
            } elsif ( $type eq 'NJ' ) {
               $BS_count = 0;
            }
            $the_input_newick = $2;
            $the_input_newick =~ s/;?\s*$//; # remove final ; and whitespace if present
         }

         #	print STDERR $the_input_newick, "\n"; exit;
         #  $next_line =~ s/;?\s*$//;	# remove final ; and whitespace if present\
         #	print STDERR "$BS_count $max_nbs  $the_input_newick  ", scalar keys %$seqid_species, " \n\n\n\n"; exit;
         if ( $BS_count <= $max_nbs ) {
            #	print STDERR "$BS_count $max_nbs  $the_input_newick  ", scalar keys %$seqid_species, " \n\n\n\n"; 
            #	print STDERR "XXX: [", ! ($the_input_newick =~ /\[species=[a-zA-Z_]+\]/), "]\n"; exit;
            if ((! ($the_input_newick =~ /\[species=[a-zA-Z_]+\]/)) and (scalar keys %$seqid_species > 0)) {
               #	print STDERR "$BS_count $the_input_newick \n\n\n\n";
               #	print "XXX\n";
               $the_input_newick = taxonify_newick( $the_input_newick, $seqid_species );
               #	print STDERR "YYY \n";
            }
	$the_input_newick =~ s/[|]/_/g; #replace pipes with underscores
            ################## CONSTRUCT TREE OBJECT #################################
            my $parser = CXGN::Phylo::Parse_newick->new( $the_input_newick, 0 );
            my $tree = $parser->parse( CXGN::Phylo::BasicTree->new() );
            $tree->impose_branch_length_minimum();
            $tree->show_newick_attribute("species");
            $tree->set_show_standard_species(0);
            $tree->get_root()->recursive_implicit_names();
            $tree->get_root()->recursive_implicit_species();

            # if species is empty, (because no [species=x] in newick), get it from the name...
            $tree->set_missing_species_from_names();

            $tree->set_show_standard_species(1);
            $tree->set_species_standardizer( CXGN::Phylo::Species_name_map->new() );

            ##################### REROOTING ##########################################
            my @leaves = $tree->get_leaves();
            my $query_node = undef;
            for my $a_leaf (@leaves) {
           #    print "X: $query_sequence_id  ", $a_leaf->get_name(), "\n";
               if ($a_leaf->get_name() eq $query_sequence_id) {
                  $query_node = $a_leaf;
                  last;
               }
            }
           # print STDERR "fam sizes: $family_size  ", scalar @leaves, "\n";
            if(($family_size != scalar @leaves) or (!defined $query_node)){
               print STDERR "Query $query_sequence_id, not found in tree, or family size discrepancy: $family_size  ", scalar @leaves, ". Skipping. \n";
               next;
            }
            my $seqid_presence = {$query_sequence_id => 1}; # this hash holds the ids in the maximal part of tree containing query and only AM hosts.
            $tree->reset_root_to_point_on_branch($query_node, 0.5*$query_node->get_branch_length());
            $tree->get_root()->recursive_implicit_names(); # is this needed?

            #   print $tree->generate_newick(), "\n";

            my $species_count = {$seqid_species->{$query_sequence_id} => 1};
            #species_and_category_counts($tree->get_root(), $seqid_species, $species_category);
            # for (keys %$species_count) {
            #    print "$query_sequence_id   $_  ", $species_count->{$_}, "\n";
            # } 
            # print "$query_sequence_id   ";
            #   for ('7_nonangiosperms', '13_nonAM_angiosperms', '35_AM_angiosperms') {
            #      my $csc = (exists $cat_spcount->{$_})? $cat_spcount->{$_} : 0;
            #                print "$csc  ";
            #             } print "\n";
            my @children = $tree->get_root()->get_children();
            
            my $next_node = $children[1];
            while (1) {
               @children = $next_node->get_children();
               my $n_children = scalar @children;

               # if($n_children > 2);
               if ($n_children == 2) {
                  my ($Lspecies_count, $Lcat_spcount, $Lid_presence) =  species_and_category_counts($children[0], $seqid_species, $species_category);
                  #     print "L:  ", join(", ", keys %$Lspecies_count), "\n";
                  my ($Rspecies_count, $Rcat_spcount, $Rid_presence) =  species_and_category_counts($children[1], $seqid_species, $species_category);
                  #     print "R:  ", join(", ", keys %$Rspecies_count), "\n";
                  my $L_AM_OK = valuez($Lcat_spcount, '13_nonAM_angiosperms') <= $max_allowed_nonAM; # number of nonAMang species found in L subtree is <= max_allowed_nonAM
                  # (!exists $Lcat_spcount->{'13_nonAM_angiosperms'} or $Lcat_spcount->{'13_nonAM_angiosperms'} <= $max_allowed_nonAM)? 1 : 0;
                  my $L_nonang_OK = valuez($Lcat_spcount, '7_nonangiosperms') <= $max_allowed_nonang;
                  my $R_AM_OK = valuez($Rcat_spcount, '13_nonAM_angiosperms') <= $max_allowed_nonAM;
                  my $R_nonang_OK = valuez($Rcat_spcount, '7_nonangiosperms') <= $max_allowed_nonang;
                  # (!exists $Rcat_spcount->{'13_nonAM_angiosperms'} or $Rcat_spcount->{'13_nonAM_angiosperms'} <= $max_allowed_nonAM)? 1 : 0;
                  #        my $L_AMangio_count = $Lcat_spcount->{'35_AM_angiosperms'};
                  #        my $R_AMangio_count = $Lcat_spcount->{'35_AM_angiosperms'};
                  my $L_OK = ($L_AM_OK and $L_nonang_OK);
                  my $R_OK = ($R_AM_OK and $R_nonang_OK);
              #    print "    [[", keys %$seqid_presence, "]] \n";
                  if ($L_OK) {
                     $seqid_presence = add_hashes($seqid_presence, $Lid_presence);
                     if ($R_OK) { # both subtrees have no nonAM's (whole tree has only AM's, nonangios)
                        ($species_count, my $cat_spcount) = species_and_category_counts($tree->get_root(), $seqid_species, $species_category);
                        $seqid_presence = add_hashes($seqid_presence, $Rid_presence);
                        last;
                     } else {   # R subtree has nonAM's
                        $species_count = add_hashes($species_count, $Lspecies_count);
                        $next_node = $children[1];
                     }
                  } else {      # L subtree has nonAM's
                     if ($R_OK) {
                        $species_count = add_hashes($species_count, $Rspecies_count);
                        $seqid_presence = add_hashes($seqid_presence, $Rid_presence);
                        $next_node = $children[0];
                     } else {   # both subtrees have nonAM's -- done
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
            }
            # for(keys %$species_count){
            #    print "$_  ", $species_count->{$_}, "\n";
            # }
            my $cat_spcount = category_count($species_count, $species_category);
            # for(keys %$cat_spcount){
            #    print "$_  ", (exists $cat_spcount->{$_})? $cat_spcount->{$_} : 0, "\n";
            # }

            my ($Amb_count, $AM_mono_count, $AM_di_count) = (valuez($cat_spcount, 'Amborella'), valuez($cat_spcount, '12_AM_monocots'), valuez($cat_spcount, '23_AM_dicots'));
            my @ids_in_qAM_subtree = keys %$seqid_presence;
            printf "$query_sequence_id    ";
            printf("%3i %3i %3i  %3i   %3i   ", 
                   $Amb_count, $AM_mono_count, $AM_di_count, 
                   sum($Amb_count, $AM_mono_count, $AM_di_count), valuez($cat_spcount, '7_nonangiosperms'));
            print "[", join(" ", @ids_in_qAM_subtree), "]", "\n";
           
            $tree->impose_branch_length_minimum(1);
            $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
         }
      } elsif ( $next_line =~ /^\s*$/ ) {
         # line has only whitespace. do nothing
      } else {
         warn "line: \n[$_] has unexpected format.";
      }
   }
   
}
close $fh_in;
# end of main

sub valuez{ # returns value of key/value pair, or 0 if pair doesn't exist for that key.
   my $href = shift;
   my $key = shift;
   return (exists $href->{$key})? $href->{$key} : 0;
}

sub category_count{
   my $sp_count = shift;
   my $sp_cat = shift;
   my $cat_spcount = {};
   for my $sp (keys %$sp_count) {
      my $cat = $sp_cat->{$sp};
      $cat_spcount->{$cat}++;
   }
   return $cat_spcount;
}

sub add_hashes{
   my $sc1 = shift;
   my $sc2 = shift;
   while ( my ($sp, $c) = each %$sc2) {
      $sc1->{$sp} += $c;
   }
   return $sc1;
}

sub species_and_category_counts{
   my $subtree_root = shift;
   my $id_species = shift;
   my $species_category = shift; # keys are species names (e.g. 'Aquilegia_coerulea'), values categories (e.g. '35_AM_angiosperms')
   my $species_count = {};
   my $category_speciescount = {};
   my %subtree_id_presence = (); # keys: ids in subtree, values: 1
   #  my $category_count = {};
   # my $node_name = $subtree_root->get_name();
   my @ids = @{$subtree_root->get_implicit_names()};
   for my $an_id (@ids) {
      $subtree_id_presence{$an_id} = 1;
#	print "id: $an_id \n";
      my $species = $id_species->{$an_id};
 #    print "id, species: $an_id  $species \n";
      my $category = $species_category->{$species};
      # print "cat: $category \n";
      if (! exists $species_count->{$species}) {
         $category_speciescount->{$category}++;
      }
      $species_count->{$species}++;
      #   $category_count->{$category}++;
   }
   # for(keys %$species_count){
   #    $category_speciescount->{$species_category->{$_}}++;
   # }
   return ($species_count, $category_speciescount, \%subtree_id_presence); # , $category_count);
}

sub reroot {
   my $tree = shift;
   my $reroot_method = shift || 'mindl';

   #    my $species_tree_filename = shift; # file name of species tree newick file. overrides default species tree
   my $species_tree = shift;
   my ( $new_root, $dist_above ) = ( undef, undef );
   if ( $reroot_method eq 'none' ) {
      return $tree;
   } elsif ( $reroot_method eq 'midpoint' ) {
      ( $new_root, $dist_above ) = $tree->find_midpoint();
   } elsif ( $reroot_method eq 'minvar' ) {
      ( $new_root, $dist_above ) = $tree->min_leaf_dist_variance_point();
   } elsif ( $reroot_method eq 'mindl' ) {

      #         my $species_tree_newick;
      # #	print STDERR "species tree filename: [$species_tree_filename]\n"; # exit;
      #         if ( defined $species_tree_filename ) {
      # #	  $species_tree_filename = abs_path($species_tree_filename); #
      # 	  print STDERR "sp. tr. filename: $species_tree_filename \n";
      #             my $species_tree_file_obj = CXGN::Phylo::File->new($species_tree_filename);
      #             $species_tree_newick = $species_tree_file_obj->get_tree_string();
      #         }
      #         else {
      #             $species_tree_newick =
      # '( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thelungiella_h[species=Thelungiella_halophila]:0.01, Thelungiella_s[species=Thelungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

      # # "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
      #         }
      #         my $species_tree_parser = CXGN::Phylo::Parse_newick->new( $species_tree_newick, 1 );
      #         my $species_tree = $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

      #         $species_tree->set_missing_species_from_names();    # get species from name if species undef
      #         $species_tree->impose_branch_length_minimum();
      #         $species_tree->collapse_tree();
      #         $species_tree->get_root()->recursive_implicit_names();
      #         $species_tree->get_root()->recursive_implicit_species();

      my $spec_bit_hash = $tree->get_species_bithash($species_tree);

      $tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);
      $species_tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);

      ( $new_root, $dist_above ) = $tree->find_mindl_node($species_tree);
   } else {
      warn "reroot option $reroot_method unknown. No rerooting performed.\n";
   }
   if ( defined $new_root ) {
      $tree->reset_root_to_point_on_branch( $new_root, $dist_above );
   } else {
      warn "In reroot. \$new_root undefined, can't reroot. Tree unchanged.\n";
   }
   return $tree;
}

sub store_gg_info {
   my $gg_filename   = shift;
   my %seqid_species = ();
   my %species_count = ();
   if ( defined $gg_filename and -f $gg_filename ) {
      open my $fh_gg, "<", "$gg_filename";
      while (<$fh_gg>) {
         my @cols = split( " ", $_ );
         my $species = shift @cols;
         $species =~ s/:$//;    # remove final colon if present.
         $species_count{$species}++;
         for (@cols) {
		my $id = $_;
	$id =~ s/[|]/_/g; # change pipes to underscores.
            $seqid_species{$id} = $species;
         }
      }
   }                # done storing gg_file info in hash %seqid_species
   return (\%seqid_species, \%species_count);
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

sub get_species_tree {
   my $species_tree_filename = shift;

   my $species_tree_newick;

   #	print STDERR "species tree filename: [$species_tree_filename]\n"; # exit;
   if ( defined $species_tree_filename ) {

      #	  $species_tree_filename = abs_path($species_tree_filename); #
      #	  print STDERR "sp. tr. filename: $species_tree_filename \n";
      my $species_tree_file_obj = CXGN::Phylo::File->new($species_tree_filename);
      $species_tree_newick = $species_tree_file_obj->get_tree_string();
   } else {
      $species_tree_newick =
        '( ((ostreococcus_t[species=Ostreococcus_tauri]:1, ostreococcus_l[species=Ostreococcus_lucimarinus]:1):1, (chlamydomonas[species=Chlamydomonas_reinhardtii]:1, volvox[species=Volvox_carteri]:1):1):1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( (norway_spruce[species=Picea_abies]:1, loblolly_pine[species=Pinus_taeda]:1):1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, (banana[species=Musa_acuminata]:1,( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( nelumbo[species=Nelumbo_nucifera]:1, ( ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, (monkey_flower[species=Mimulus_guttatus]:1, bladderwort[species=Utricularia_gibba]:1):1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, beet[species=Beta_vulgaris]:1):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( cleome[species=Tarenaya_hassleriana]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1):1, ( ( ( lupine[species=Lupinus_angustifolius]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1)';


      #'( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thelungiella_h[species=Thelungiella_halophila]:0.01, Thelungiella_s[species=Thelungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

      # "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
   }
   my $species_tree_parser = CXGN::Phylo::Parse_newick->new( $species_tree_newick, 0 );
   my $species_tree = $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

   $species_tree->set_missing_species_from_names(); # get species from name if species undef
   $species_tree->impose_branch_length_minimum();
   $species_tree->collapse_tree();
   $species_tree->get_root()->recursive_implicit_names();
   $species_tree->get_root()->recursive_implicit_species();
   return $species_tree;
}
