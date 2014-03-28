#!/usr/bin/perl -w
use strict;


use lib '/home/tomfy/cxgn/cxgn-corelibs/lib/';

use CXGN::Phylo::Overlap;
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;


my $default_newick = '((((((Bradi1g50820.1[species=Brachypodium_distachyon]:0.01,Bradi1g48660.1[species=Brachypodium_distachyon]:0.01):0.0001,(Bradi4g07840.1[species=Brachypodium_distachyon]:0.01,Bradi2g27720.1[species=Brachypodium_distachyon]:0.01):0.0001):0.0001,((Bradi2g18410.1[species=Brachypodium_distachyon]:0.01,Bradi2g24080.1[species=Brachypodium_distachyon]:0.01):0.0001,(Medtr8g103245.1[species=Medicago_truncatula]:0.01,Medtr8g092820.1[species=Medicago_truncatula]:0.01):0.0001):0.0001):0.0001,(((Medtr8g092720.1[species=Medicago_truncatula]:0.01,Medtr8g063500.1[species=Medicago_truncatula]:0.01):0.0001,(Medtr7g059070.1[species=Medicago_truncatula]:0.01,Medtr5g029820.1[species=Medicago_truncatula]:0.01):0.0001):0.0001,((Medtr5g029770.1[species=Medicago_truncatula]:0.01,Medtr4g088150.1[species=Medicago_truncatula]:0.01):0.0001,(Medtr4g065990.1[species=Medicago_truncatula]:0.01,Medtr2g035230.1[species=Medicago_truncatula]:0.01):0.0001):0.0001):0.0001):0.0001,((((POPTR_0005s25530.1[species=Populus_trichocarpa]:0.01,Glyma18g50120.1[species=Glycine_max]:0.01):0.0001,(Glyma15g34870.1[species=Glycine_max]:0.01,Glyma15g03710.1[species=Glycine_max]:0.01):0.0001):0.0001,((Glyma13g41710.1[species=Glycine_max]:0.01,Glyma13g24180.1[species=Glycine_max]:0.01):0.0001,(Glyma12g05910.1[species=Glycine_max]:0.01,Glyma11g13940.1[species=Glycine_max]:0.01):0.0001):0.0001):0.0001,(((Glyma08g26880.1[species=Glycine_max]:0.01,Glyma08g24380.1[species=Glycine_max]:0.01):0.0001,(Glyma05g32200.1[species=Glycine_max]:0.01,Cucsa.365170.1[species=Cucumis_sativus]:0.01):0.0001):0.0001,((GRMZM2G387076_P03[species=Zea_mays]:0.01,GRMZM2G451254_P01[species=Zea_mays]:0.01):0.0001,((Glyma05g32210.1[species=Glycine_max]:0.01,Glyma07g32360.1[species=Glycine_max]:0.01):0.0001,(Medtr7g013610.1[species=Medicago_truncatula]:0.01,Bradi3g45290.1[species=Brachypodium_distachyon]:0.01):0.0001):0.0001):0.0001):0.0001):0.0001):0.001905,((((((Sb06g016850.1[species=Sorghum_bicolor]:0.01,Sb06g016330.1[species=Sorghum_bicolor]:0.01):0.0001,(Sb04g022160.1[species=Sorghum_bicolor]:0.01,Sb03g005550.1[species=Sorghum_bicolor]:0.01):0.0001):0.0001,((LOC_Os11g05730.1[species=Oryza_sativa]:0.01,LOC_Os06g06510.1[species=Oryza_sativa]:0.01):0.0001,(LOC_Os06g06460.1[species=Oryza_sativa]:0.01,LOC_Os05g36280.1[species=Oryza_sativa]:0.01):0.0001):0.0001):0.0001,(((LOC_Os01g64640.1[species=Oryza_sativa]:0.01,Aly471009[species=Arabidopsis_lyrata]:0.01):0.0001,(Aly892824[species=Arabidopsis_lyrata]:0.01,Aly908994[species=Arabidopsis_lyrata]:0.01):0.0001):0.0001,((Aly484525[species=Arabidopsis_lyrata]:0.01,Aly487869[species=Arabidopsis_lyrata]:0.01):0.0001,(Thhalv10014985m[species=Thellungiella_halophila]:0.01,Thhalv10014986m[species=Thellungiella_halophila]:0.01):0.0001):0.0001):0.0001):0.0001,((((Thhalv10005111m[species=Thellungiella_halophila]:0.01,Thhalv10005114m[species=Thellungiella_halophila]:0.01):0.0001,(Thhalv10029057m[species=Thellungiella_halophila]:0.01,Carubv10003405m[species=Capsella_rubella]:0.01):0.0001):0.0001,((Carubv10002544m[species=Capsella_rubella]:0.01,Carubv10011755m[species=Capsella_rubella]:0.01):0.0001,(Carubv10018763m[species=Capsella_rubella]:0.01,evm.model.supercontig_84.95[species=Carica_papaya]:0.01):0.0001):0.0001):0.0001,(((evm.model.supercontig_84.93[species=Carica_papaya]:0.01,evm.model.supercontig_84.92[species=Carica_papaya]:0.01):0.0001,(POPTR_0002s03030.1[species=Populus_trichocarpa]:0.01,POPTR_0001s05470.1[species=Populus_trichocarpa]:0.01):0.0001):0.0001,((POPTR_0014s09260.1[species=Populus_trichocarpa]:0.01,POPTR_0003s22240.1[species=Populus_trichocarpa]:0.01):0.0001,(POPTR_0003s22120.1[species=Populus_trichocarpa]:0.01,Cucsa.143240.1[species=Cucumis_sativus]:0.01):0.0001):0.0001):0.0001):0.0001):0.0001,(((((Cucsa.065360.1[species=Cucumis_sativus]:0.01,Cucsa.050120.1[species=Cucumis_sativus]:0.01):0.0001,(Cucsa.050110.1[species=Cucumis_sativus]:0.01,PGSC0003DMP400055010[species=Solanum_tuberosum]:0.01):0.0001):0.0001,((PGSC0003DMP400011232[species=Solanum_tuberosum]:0.01,PGSC0003DMP400055218[species=Solanum_tuberosum]:0.01):0.0001,(PGSC0003DMP400039599[species=Solanum_tuberosum]:0.01,PGSC0003DMP400011206[species=Solanum_tuberosum]:0.01):0.0001):0.0001):0.0001,(((PGSC0003DMP400043844[species=Solanum_tuberosum]:0.01,Solyc10g008910.1.1[species=Solanum_lycopersicum]:0.01):0.0001,(Solyc02g077480.1.1[species=Solanum_lycopersicum]:0.01,Solyc01g086820.2.1[species=Solanum_lycopersicum]:0.01):0.0001):0.0001,((Solyc01g080600.2.1[species=Solanum_lycopersicum]:0.01,Solyc01g074000.2.1[species=Solanum_lycopersicum]:0.01):0.0001,(Solyc01g073970.2.1[species=Solanum_lycopersicum]:0.01,GRMZM2G376957_P01[species=Zea_mays]:0.01):0.0001):0.0001):0.0001):0.0001,((((GRMZM2G401581_P01[species=Zea_mays]:0.01,GRMZM2G130079_P02[species=Zea_mays]:0.01):0.0001,(GRMZM2G130079_P01[species=Zea_mays]:0.01,GRMZM5G864735_P01[species=Zea_mays]:0.01):0.0001):0.0001,((GRMZM2G418258_P01[species=Zea_mays]:0.01,GRMZM2G179005_P03[species=Zea_mays]:0.01):0.0001,(GRMZM2G179005_P02[species=Zea_mays]:0.01,GRMZM2G355773_P01[species=Zea_mays]:0.01):0.0001):0.0001):0.0001,(((AT5G65360.1[species=Arabidopsis_thaliana]:0.01,AT5G10400.1[species=Arabidopsis_thaliana]:0.01):0.0001,(AT5G10390.1[species=Arabidopsis_thaliana]:0.01,AT3G27360.1[species=Arabidopsis_thaliana]:0.01):0.0001):0.0001,((AT1G09200.1[species=Arabidopsis_thaliana]:0.01,29900.m001544[species=Ricinus_communis]:0.01):0.0001,((Carubv10021660m[species=Capsella_rubella]:0.01,Sb10g004110.1[species=Sorghum_bicolor]:0.01):0.0001,(Sb10g004100.1[species=Sorghum_bicolor]:0.01,Sb09g021650.1[species=Sorghum_bicolor]:0.01):0.0001):0.0001):0.0001):0.0001):0.0001):0.0001):0.005715)';

my $newick = shift || $default_newick;;
my  $species_tree_newick =
      '( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';
  my $species_tree_parser =
    CXGN::Phylo::Parse_newick->new( $species_tree_newick, 0 );
  my $species_tree =
    $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

  $species_tree->set_missing_species_from_names()
    ;			      # get species from name if species undef
  $species_tree->impose_branch_length_minimum();
  $species_tree->collapse_tree();
  $species_tree->get_root()->recursive_implicit_names();
  $species_tree->get_root()->recursive_implicit_species();

print STDERR "ZZZZZ\n";

my $rrnewick = newick2rerootednewick($newick, 'mindl', $species_tree);

print STDERR "RRR: $rrnewick \n";

sub newick2rerootednewick
  { # take a taxonified newick expression, make tree, reroot is by mindl wrt $species_tree

    # return newick expression
    my $taxonified_newick = shift;
$taxonified_newick = put_min_branchlengths_into_newick($taxonified_newick, 0.01);
print STDERR "ZZZ: $taxonified_newick \n";
    my $reroot_method     = shift;
    my $species_tree      = shift;
    my $parser = CXGN::Phylo::Parse_newick->new( $taxonified_newick, 0 );
   print STDERR "after constructing parser \n";
    my $tree   = $parser->parse( CXGN::Phylo::BasicTree->new() );
    $tree->impose_branch_length_minimum();
    $tree->show_newick_attribute("species");
 print STDERR "AAAA:  ", $tree->generate_newick(), "\n";
    $tree->set_show_standard_species(0);
    $tree->get_root()->recursive_implicit_names();
    $tree->get_root()->recursive_implicit_species();
 print STDERR "AAAA:  ", $tree->generate_newick(), "\n";
    # if species is empty, (because no [species=x] in newick), get it from the name...
    #$tree->set_missing_species_from_names();

    $tree->set_show_standard_species(1);
    $tree->set_species_standardizer( CXGN::Phylo::Species_name_map->new() );
 print STDERR "AAAA:  ", $tree->generate_newick(), "\n";
    ### make the tree binary
    my @new_root_point;
    { # if root has > 2 children, reroot to pt on a root-child branch, so new root has 2 children.
      my @root_children = $tree->get_root()->get_children();
      if ( scalar @root_children > 2 ) {
	@new_root_point = (
			   $root_children[0], 0.9 * $root_children[0]->get_branch_length()
			  );
	$tree->reset_root_to_point_on_branch(@new_root_point);
      }
    }

    # binarify every non-binary node. At present doesn't attempt to choose in a smart way
    # among the various possible resolutions
 print STDERR "AAAA:  ", $tree->generate_newick(), "\n";
    $tree->make_binary( $tree->get_root() ); # urec requires binary tree.
print STDERR "Before reroot.\n";
print STDERR "newick: \n", $tree->generate_newick(), "\n\n";
    $tree = reroot( $tree, $reroot_method, $species_tree );
    my $rerooted_tree_newick = $tree->generate_newick();

    $tree->decircularize();    # so can be properly garbage collected!
    return $rerooted_tree_newick;
  }

sub put_min_branchlengths_into_newick{
  my $newick = shift;
  my $min_bl = shift || 0.00001;
$newick =~ s/:0[.]0+([,)])/:$min_bl$1/g;
return $newick;
}

sub reroot {
  my $tree             = shift;
  my $reroot_method    = shift || 'mindl';
  my $species_tree     = shift;
  my $species_tree_arg = shift; # file name of species tree newick file. overrides default species tree

  #  print "in reroot. reroot method $reroot_method \n";
  #print $tree->generate_newick(), "\n";
  my ( $new_root, $dist_above ) = ( undef, undef );
  if ( $reroot_method eq 'none' ) {
    return $tree;
  } elsif ( $reroot_method eq 'midpoint' ) {
    ( $new_root, $dist_above ) = $tree->find_midpoint();
  } elsif ( $reroot_method eq 'minvar' ) {
    ( $new_root, $dist_above ) = $tree->min_leaf_dist_variance_point();
  } elsif ( $reroot_method eq 'mindl' ) {

    #   if(0){
    #   my $species_tree_newick;
    #   if ( defined $species_tree_arg ) {
    #     my $species_tree_file = CXGN::Phylo::File->new($species_tree_arg);
    #     find_cycle($species_tree_file);
    #     $species_tree_newick = $species_tree_file->get_tree_string();
    #   } else {
    #     $species_tree_newick =
    # 	'( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

    #     # "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
    #   }
    #   my $species_tree_parser = CXGN::Phylo::Parse_newick->new( $species_tree_newick, 1 );
    #   my $species_tree = $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

    #   $species_tree->set_missing_species_from_names(); # get species from name if species undef
    #   $species_tree->impose_branch_length_minimum();
    #   $species_tree->collapse_tree();
    #   $species_tree->get_root()->recursive_implicit_names();
    #   $species_tree->get_root()->recursive_implicit_species();
    # }

    my $spec_bit_hash = $tree->get_species_bithash($species_tree);

    $tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);
    $species_tree->get_root()
      ->recursive_set_implicit_species_bits($spec_bit_hash);

    ( $new_root, $dist_above ) = $tree->find_mindl_node($species_tree);

    #  print "after find_mindl. new root node name: ", $new_root->get_name(), " $dist_above\n";
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
