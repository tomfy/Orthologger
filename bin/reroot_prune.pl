#!/usr/bin/perl -w
use strict;
use warnings;
no lib '/home/tomfy/bin';
no lib '/home/tomfy/Orthologger/bin';
no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {
    $bindir = dirname( abs_path(__FILE__) );    # this has to go in Begin block so happens at compile time
    $libdir = $bindir . '/../lib';
}
use lib $libdir;

use Getopt::Long;

use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

my $default_gg_file_path           = "$bindir/../tst/21species.gg";
my $default_species_tree_file_path = "$bindir/../species_tree_plant_52.newick";

# Defaults
my $input_newick_filename;
my $sequence_id;
my $gg_filename         = $default_gg_file_path;
my $do_ortholog_support = 0;                       # shift || 0; # otherwise just get tree with FastTree
my $do_taxonify     = 1;
my $reroot_method   = 'mindl';
my $prune_threshold = 3;

my $clade_outer_edge = 'monocots';    # default is to expand clade until contains $prune_threshold  monocot sequences
my $species_tree_newick_file = $default_species_tree_file_path;
if ( !-f $species_tree_newick_file ) {
    die "$species_tree_newick_file is not a regular file. Will use default species tree.\n";
    $species_tree_newick_file = undef;
}
my $others = 'Arabidopsis_thaliana';

# Process long cl options
GetOptions(
    'input_newick=s' => \$input_newick_filename,
    'sequence_id=s'  => \$sequence_id,
    'ggfile=s'       => \$gg_filename
    , # defines species-sequence association. 1 line per species: species name followed by whitespace-separated sequence ids.
    'reroot=s' => \$reroot_method,    # selects rerooting method. options: none, midpoint, minvar, mindl.
    'prune_threshold=i' =>
      \$prune_threshold,    # prune to minimal tree containing id of interest, and $prune_threshold non-dicot sequences.
    'speciestreefile=s' => \$species_tree_newick_file
    ,                       # to override built-in species tree with 52 species (see sub reroot below for this tree).
    'clade_outer_edge=s' => \$clade_outer_edge,
    'other_species=s'    => \$others,
);

my $the_input_newick = `cat $input_newick_filename`;
my $taxonified_newick = taxonify_newick($the_input_newick, $gg_filename);
#print "[$the_input_newick]\n";
################## CONSTRUCT TREE OBJECT #################################
my $parser = CXGN::Phylo::Parse_newick->new( $taxonified_newick, 1 );
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

#print "Tree object constructed.\n";
##################### REROOTING ##########################################
#    my $opt_r = 'midpoint';
#    my $custom_species_tree_newick_file = undef; # needed only for mindl rerooting. and reroot has default tree for this.
$tree = reroot( $tree, $reroot_method, $species_tree_newick_file );

my $whole_tree_newick = $tree->generate_newick();

#  print "whole rerooted tree: $whole_tree_newick \n";
my $newick_filename = $input_newick_filename;

# $newick_filename =~ s/[.]fasta$//;
# $newick_filename =~ s/_align//;
# $newick_filename =~ s/_tax//;
$newick_filename =~ s/[.]newick$//;
my $pruned_newick_filename = $newick_filename . "_rerooted_pruned.newick";
$newick_filename .= "_rerooted.newick";
open my $fh_whole_newick, ">$newick_filename";
print $fh_whole_newick "$whole_tree_newick\n";
close $fh_whole_newick;

############################# PRUNING ####################################
my $non_angiosperm_species = {
    'Chlamydomonas_reinhardtii'  => 1,
    'Physcomitrella_patens'      => 1,
    'Selaginella_moellendorffii' => 1,
    'Pinus_taeda'                => 1,
};
my $amborella = { 'Amborella_trichopoda' => 1 };

my $monocot_species = {
    'Phoenix_dactylifera' => 1,
    'Setaria_italica'     => 1,
    'Triticum_aestivum'   => 1,

    'Hordeam_vulgare'         => 1,
    'Zea_mays'                => 1,
    'Brachypodium_distachyon' => 1,
    'Sorghum_bicolor'         => 1,
    'Oryza_sativa'            => 1
};
my $clade_outer_edge_species;

if ( $clade_outer_edge eq 'monocots' ) {
    $clade_outer_edge_species = $monocot_species;
}
else {
    my $user_outer_edge_species = {};
    my @coe_species = split( "[, ]+", $clade_outer_edge ); # separated by comma, whitespace or both
    for (@coe_species) {
        $user_outer_edge_species->{$_} = 1;
    }
    $clade_outer_edge_species = $user_outer_edge_species;
}

print "Prune to the minimal clade containing sequence: $sequence_id, \n";
print "and at least $prune_threshold sequences whose species are in the following list: \n";
print keys_string( $clade_outer_edge_species, 22, 4 );

my $other_species;
my @oth_species = split( "[, ]+", $others );
for (@oth_species) {
    $other_species->{$_} = 1;
}
print "Other species to count in clade and out of clade: \n";
print keys_string( $other_species, 22, 4 );

my ( $pruned_tree_root, $in_other, $out_other ) =
  $tree->min_clade( $sequence_id, $prune_threshold, $clade_outer_edge_species, $other_species );
print "Other species counts; in clade: $in_other, out of clade: $out_other \n";
my $pruned_tree_newick = $pruned_tree_root->recursive_generate_newick();

#  print "pruned tree: $pruned_tree_newick \n";
#    my $pruned_newick_filename = $align_filename . "_pruned.newick";
open my $fh_pruned_newick, ">$pruned_newick_filename";
print $fh_pruned_newick "$pruned_tree_newick\n";
close $fh_pruned_newick;

sub reroot {
    my $tree             = shift;
    my $reroot_method    = shift || 'mindl';
    my $species_tree_arg = shift;                # file name of species tree newick file. overrides default species tree

    #  print "in reroot. reroot method $reroot_method \n";
    #print $tree->generate_newick(), "\n";
    my ( $new_root, $dist_above ) = ( undef, undef );
    if ( $reroot_method eq 'none' ) {
        return $tree;
    }
    elsif ( $reroot_method eq 'midpoint' ) {
        ( $new_root, $dist_above ) = $tree->find_midpoint();
    }
    elsif ( $reroot_method eq 'minvar' ) {
        ( $new_root, $dist_above ) = $tree->min_leaf_dist_variance_point();
    }
    elsif ( $reroot_method eq 'mindl' ) {

        my $species_tree_newick;
        if ( defined $species_tree_arg ) {
            my $species_tree_file = CXGN::Phylo::File->new($species_tree_arg);
            $species_tree_newick = $species_tree_file->get_tree_string();
        }
        else {
            $species_tree_newick =
'( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thelungiella_h[species=Thelungiella_halophila]:0.01, Thelungiella_s[species=Thelungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

# "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
        }
        my $species_tree_parser = CXGN::Phylo::Parse_newick->new( $species_tree_newick, 1 );
        my $species_tree = $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

        $species_tree->set_missing_species_from_names();    # get species from name if species undef
        $species_tree->impose_branch_length_minimum();
        $species_tree->collapse_tree();
        $species_tree->get_root()->recursive_implicit_names();
        $species_tree->get_root()->recursive_implicit_species();

        my $spec_bit_hash = $tree->get_species_bithash($species_tree);

        $tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);
        $species_tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);

        ( $new_root, $dist_above ) = $tree->find_mindl_node($species_tree);

        #  print "after find_mindl. new root node name: ", $new_root->get_name(), " $dist_above\n";
    }
    else {
        warn "reroot option $reroot_method unknown. No rerooting performed.\n";
    }
    if ( defined $new_root ) {
        $tree->reset_root_to_point_on_branch( $new_root, $dist_above );
    }
    else {
        warn "In reroot. \$new_root undefined, can't reroot. Tree unchanged.\n";
    }
    return $tree;
}

sub keys_string {
    my $href      = shift;
    my $col_width = shift || 24;
    my $cols      = shift || 4;
    my $line_initial_whitespace = shift || '  ';
    my $n_keys    = scalar keys %$href;
    my $count     = 0;
    my $string    = $line_initial_whitespace;
    my $format    = "%-$col_width" . 's  ';
    while ( my ( $k, $v ) = each %$href ) {

        if ( defined $v ) {
            $string .= sprintf( "$format", $k );
            $count++;

            #    print "count, cols, nkeys: $count, $cols, $n_keys \n";
            $string .= sprintf("\n$line_initial_whitespace") if ( $count % $cols == 0 or $count == $n_keys );
        }
        else {
            $n_keys--;
        }
    }
    $string =~ s/$line_initial_whitespace$//;
    return $string;
}


sub taxonify_newick{
  my $newick = shift;
  my $gg_filename = shift;
  my %seqid_species = ();
  if (defined $gg_filename and -f $gg_filename) {
    open my $fh_gg, "<", "$gg_filename";
    while (<$fh_gg>) {
      my @cols = split(" ", $_);
      my $species = shift @cols;
      $species =~ s/:$//;	# remove final colon if present.
      for (@cols) {
	$seqid_species{$_} = $species;
      }
    }
  }	      # done storing gg_file info in hash %seqid_species
$newick =~ s/\s+$//; # remove final whitespace
  my $new_newick = $newick;
  $new_newick =~ s/ ([\(,]) \s* ([^,\):]+) \s* ([,\):])/$1$2\n$3/xg; # add newline after each leaf
  my @leaf_lines = split("\n", $new_newick);
  my $last_bit = pop @leaf_lines;
  for (@leaf_lines) {
    if (/\[species=.*\]/) { # species info already there for this leaf - leave it alone
    } else {
#      print "XXX: $_ \n";
      / [\(,] \s* ([^\(,\):]+) \s* $/x; 
      my $seq_id = $1;
#      print "YYY: $seq_id \n";
      if (exists $seqid_species{$seq_id}) {
	my $species = $seqid_species{$seq_id};
	$seq_id .= '[species=' . $species . ']';
      } else {
	warn "sequence $seq_id; no corresponding species found.\n";
      }
      s/ ([\(,]) \s* ([^\(,\):]+) \s* $/$1$seq_id/x;  
    }
  }
  $new_newick = join('', @leaf_lines) . $last_bit;
  return $new_newick;
}
