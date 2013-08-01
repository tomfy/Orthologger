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

my $default_gg_file_path           = undef; #"$bindir/../tst/21species.gg";
my $default_species_tree_file_path = undef; #"$bindir/../species_tree_plant_52.newick";

# Defaults:

my $input_newicks_filename;
my $gg_filename   = $default_gg_file_path;
my $reroot_method = undef; # default is do not reroot ( mindl rerooting is default for fttree )

my $clade_specifiers = '7dicots,6 : 4monocots,3 : Selaginella_moellendorffii,1';

# This would mean we will look for a clade with 6 out of 7 dicots, for one with 3 monocots, and for one with Selaginella

my $predefined_taxon_groups = {
			       '4nonangiosperms' => {
						     'Chlamydomonas_reinhardtii'  => 1,
						     'Physcomitrella_patens'      => 1,
						     'Selaginella_moellendorffii' => 1,
						     'Pinus_taeda'                => 1,
						    },
			       '8monocots' => {
					       'Phoenix_dactylifera'     => 1, # data palm
					       'Setaria_italica'         => 1, # foxtail millet
					       'Triticum_aestivum'       => 1, # wheat
					       'Hordeam_vulgare'         => 1, # barley
					       'Zea_mays'                => 1, # maize
					       'Brachypodium_distachyon' => 1,
					       'Sorghum_bicolor'         => 1,
					       'Oryza_sativa'            => 1 # rice
					      },
			       '4monocots' => {
					       'Zea_mays'                => 1, # maize
					       'Brachypodium_distachyon' => 1,
					       'Sorghum_bicolor'         => 1,
					       'Oryza_sativa'            => 1 # rice
					      },

			       '7dicots' => {
					     'Solanum_lycopersicum' => 1, # tomato
					     'Solanum_tuberosum'    => 1, # potato
					     'Vitis_vinifera'       => 1, # grape
					     'Glycine_max'          => 1, # soy
					     'Populus_trichocarpa'  => 1, # poplar
					     'Ricinus_communis'     => 1, # castor
					     'Cucumis_sativus'      => 1 # cucumber
					    },

			       '5brassicas' => {
						Brassica_rapa           => 1, # turnip
						Arabidopsis_thaliana    => 1,
						Arabidopsis_lyrata      => 1,
						Thellungiella_halophila => 1,
						Capsella_rubella        => 1
					       }
			      };

# For each clade, count the number of these species present: (
my $disallowed_species = '5brassicas';
my $max_nbs = 1000000; # big number - do all bs replicates by default.
my $species_tree_newick_file = $default_species_tree_file_path;	# if undefined or


# Process long cl options
GetOptions(
	   'input_newicks=s' => \$input_newicks_filename,
	   'ggfile=s'        => \$gg_filename, # defines species-sequence association.
	   # gg file has 1 line per species: species name followed by whitespace-separated sequence ids.
	   'reroot=s' => \$reroot_method, # selects rerooting method. options: none, midpoint, minvar, mindl.
	   'speciestreefile=s'    => \$species_tree_newick_file, # to override built-in species tree with 52 species
	   # (see sub reroot below to see this default species tree).
	   'clade_specifiers=s'   => \$clade_specifiers,
	   'disallowed_species=s' => \$disallowed_species,
	   'maxnbs=i' => \$max_nbs,
	  );

if ( defined $species_tree_newick_file and !-f $species_tree_newick_file ) {
  die "$species_tree_newick_file is not a regular file. Will use default species tree.\n";
  $species_tree_newick_file = undef;
}
my $species_tree_obj = get_species_tree($species_tree_newick_file);
my @disallowed_species = ();
if ( exists $predefined_taxon_groups->{$disallowed_species} ) {
  @disallowed_species = keys %{ $predefined_taxon_groups->{$disallowed_species} };
} else {
  @disallowed_species = split( /\s*,\s*/, $disallowed_species );
}
###################### print clade specifiers, disallowed species. ######################
print "# $clade_specifiers \n#\n";
$clade_specifiers =~ s/\s+//g;
my @clade_specs = split( ":", $clade_specifiers );

my @clade_spec_objs = ();
for (@clade_specs) {
  push @clade_spec_objs, CXGN::Phylo::CladeSpecifier->new( $_, $predefined_taxon_groups );
}
while ( my ( $i, $cso ) = each @clade_spec_objs ) {
  print "# Clade ", $i + 1, " specs: \n", "# ", $cso->as_string(), "#\n";
}
print "# Disallowed species: ", join( ", ", @disallowed_species ), "\n#\n";
##########################################################################################

if ( !defined $input_newicks_filename or !-f $input_newicks_filename ) {
  die "No input newicks file specified, or file [$input_newicks_filename] does not exist.";
}
open my $fh_in, "<", "$input_newicks_filename"
  || die "Couldnt open $input_newicks_filename for reading.\n";
my $seqid_species = store_gg_info($gg_filename);
my $BS_count = 0;
while (<$fh_in>) {
  if (/^Id\s+(\S+)/) {
    my $sequence_id = $1;
    my $next_line   = <$fh_in>;
    my $type = undef;
#print STDERR "NEXT LINE: [", $next_line, "]\n";
    if ( $next_line =~ /^((\S+)\s+)?\s*\(/ ) { # newick string on this line
      my $the_input_newick;
      if ($next_line =~ /^\s*\(/ ) {
	$the_input_newick = $next_line;
	$type = 'ML';
#	print STDERR "A: \n", "$next_line \n";
      } elsif ($next_line =~ /^\s*(\S+)\s+(\(.*)$/) {
	$type = $1;
	if($type eq 'NJ_BS'){
	  $BS_count++;
	}elsif($type eq 'NJ'){
	  $BS_count = 0;
	}
	$the_input_newick = $2;
#	print STDERR "B: $1  \n", "B: $the_input_newick \n";
	$the_input_newick =~  s/;?\s*$//;	# remove final ; and whitespace if present
      }
#	print STDERR $the_input_newick, "\n"; exit;
    #  $next_line =~ s/;?\s*$//;	# remove final ; and whitespace if present\
           
      if($BS_count <= $max_nbs){
      if ( defined $gg_filename and -f $gg_filename ) {
	$the_input_newick = taxonify_newick( $the_input_newick, $seqid_species );
      }
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
      if ( defined $reroot_method ) {
	$tree = reroot( $tree, $reroot_method, $species_tree_obj);
      }

      ############################# FINDING SPECIFIED CLADES  ####################################
      my $cladeinfo = $tree->find_clades( $sequence_id, \@clade_spec_objs, \@disallowed_species );
      printf( "%15s %10s  ", $sequence_id, $type );
      for (@clade_spec_objs) {
	my $info = $cladeinfo->{$_};
	printf( "   %3i %3i %3i", $info->{nodes_up}, $info->{n_same}, $info->{n_disallowed_species} );
      }
      print "\n";

      $tree->impose_branch_length_minimum(1);
      $tree->decircularize(); # done with tree - decircularize so can be garbage collected.
    }
    } elsif ( $next_line =~ /^\s*$/ ) {
#      print STDERR "{{{]]]\n";
      # line has only whitespace. do nothing
    } else {
      warn "line $_ has unexpected format.";
    }
  }

}
close $fh_in;

# end of main

sub reroot {
  my $tree             = shift;
  my $reroot_method    = shift || 'mindl';
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
  if ( defined $gg_filename and -f $gg_filename ) {
    open my $fh_gg, "<", "$gg_filename";
    while (<$fh_gg>) {
      my @cols = split( " ", $_ );
      my $species = shift @cols;
      $species =~ s/:$//;	# remove final colon if present.
      for (@cols) {
	$seqid_species{$_} = $species;
      }
    }
  }		    # done storing gg_file info in hash %seqid_species
  return \%seqid_species;
}

sub taxonify_newick {
  my $newick        = shift;
  my %seqid_species = %{ shift @_ };
  $newick =~ s/\s+$//;		# remove final whitespace
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

sub get_species_tree{
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
      '( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thelungiella_h[species=Thelungiella_halophila]:0.01, Thelungiella_s[species=Thelungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

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
