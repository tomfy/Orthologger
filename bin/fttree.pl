#!/usr/bin/perl -w
use strict;

# use Devel::Cycle;

no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {
  $bindir =
    dirname( abs_path(__FILE__) ); # this has to go in Begin block so happens at compile time
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;

use CXGN::Phylo::Overlap;
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

# read in a file with alignments for each of many families
# The format is for each family:
# 1) 1 line of info about family as whole
# 2) fasta for all aligned sequences in family (this is omitted in cases of families which do not
#    meet certain criteria; at present by default if they don't include sequences from >= 3 monocot species
# 3) one blank line
# for example: (this just shows the sequence (with gaps) for the first two of the 
# Id Medtr1g004990.1 family. fam_size: 201 Arabidopsis_lyrata,Arabidopsis_thaliana,Brachypodium_distachyon,Brassica_rapa,Capsella_
# rubella,Carica_papaya,Chlamydomonas_reinhardtii,Cucumis_sativus,Glycine_max,Medicago_truncatula,Oryza_sativa,Physcomitrella_pate
# ns,Populus_trichocarpa,Ricinus_communis,Selaginella_moellendorffii,Solanum_lycopersicum,Solanum_tuberosum,Sorghum_bicolor,Thellu
# ngiella_halophila,Vitis_vinifera,Zea_mays  1 
# >Medtr1g004990.1 
# -------------------------MEKVVGGKYRIGR-KIGSGSFGEIYIGAHVVTSEL
# VAIKKEKKKTQQPQLLYEAKLYNILKGGSGIPRMKWFGTDGDYNVLVLELMGPSLDDLLY
# YCSGKFSLKSVLMLADQMLTRIEYLHSKGLLHRDIKPDNFLMGLGKKANQ--VCMIDFGL
# SKGYRDPISYKHIPYRENKNLTGTARYASSNTHKGIEQSRRDDLESLGYVLLYFLRGSLP
# WQGLQAATRMQKYEKICETKLNTPIEVLCKSCPVEFASYFHYCRSLTFDQRPDYGYLKRL
# FRELFTSKGY----------------------------AADYLYDWTILKYQE---IQ--
# QIKEQ---NQ------------SIAPVAVPTSLEPGDVDEHRE-----------------
# --------------YNDCTQNVVPKPK--------------------IYTDRPRVCMKLR
# VANVDNL--------DDEIQTDKQKVNTDLPISTVMPTED-----VPKPETTVETSNPND
# ----------------VLGSKCGASDDLVPSIRRVSSIN---------------------
# --
# >Glyma17g28670.1 
# -------------------------MERVLGGKFKVGK-KIGSGSFGEIHIGAHIETSEI
# VAIKMENRKTNQPQLQFEAKLYSTLQGGSGIPRMKWCGTDGDSNVLVIELLGPSLEDLFF
# FCGNKFSLKTVLMLADQLLTRIEYLHSKGFLHRDIKPDNFLMGLGKKANQ--VYMIDFGL
# AKEYRDPFTNKHIPYRENKGLTGTARYASYNAHSGIEQSRRDDLESLGYVLMYFLRGSLP
# WQGLQAVTKRQKYDKICKKKLSTPIEILCKSYPVEFASYFHYCRSLTFDQRPDYGLLKRL
# FRNLFTRAGY----------------------------DSDYLFDWTILKYQQ---MQ--
# QEKTQ---SQ--------------------PP----------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# --
#
# Id Medtr...
my $gg_filename = shift; 
my $n_taxa = shift || 21;
my $species_tree_newick_file = shift || undef;
my $reroot_method = shift || 'mindl';
my $min_overlap_length = 40;
my $min_taxa = 4;
my $nongap_fraction = 0.8;
my $state = 'idline';		# other possible values: fasta
my ($qid, $fam_size, $taxa, $idline, $fasta, $do);
my $support_string = '-nosupport'; # set to '' to get branch support calculation.

my $gg_hashref = store_gg_info($gg_filename);
# print STDERR "number of keys of gg_hash: ", scalar keys %$gg_hashref, "\n";

my $species_tree = get_species_tree($species_tree_newick_file);

while (<>) {
  if ($state eq 'idline') {
    my @cols = split(" ", $_);
    ($qid, $fam_size, $taxa, $do) = @cols[1,4,5,6];
    print STDERR "$qid, $fam_size, [$do]\n";
    # my @species = split(",", $taxa);
    # if(scalar @species < $min_taxa){ # need to have at least 4 taxa (Medtr + 3 monocots)
    #   $do = 0;
    # }elsif(scalar @species >= $n_taxa-1){ # if have at least 20 out of 21, guarantees 3 monocots
    #   $do = 1;
    # }else{ # need to look at taxa in more detail:
    #   my ($monocot_count, $selaginella_present) = check_taxon_list($taxa);
    #   $do = ($monocot_count >= 3)? 1 : 0;
    # }
    $idline = $_;
    $state = 'fasta';
    $fasta = '';
  } elsif ($state eq 'fasta') {
#    $fasta .= $_;
    if (/^\s*$/) { # blank line after sequence -> process the sequence.
      chomp $idline;
      my $string_to_print = "$idline   ";
      $do = ($fasta ne '');
      if ($do) {
#	print STDERR "[[[$fasta]]]\n";
	print STDERR "before overlap \n";
	my $overlap_obj = CXGN::Phylo::Overlap->new( $fasta, $nongap_fraction, 1234567 );
	my $overlap_fasta_string = $overlap_obj->overlap_fasta_string('');
	my $overlap_length       = $overlap_obj->get_overlap_length();
	print STDERR "after overlap. overlap length: $overlap_length.\n";
#print STDERR "[[[$overlap_fasta_string]]]\n";
	$string_to_print .= "$overlap_length ";
	print "$string_to_print \n";
	if($overlap_length >= $min_overlap_length){
	  my $tmpfilename = "PID" . $$ . "_tmpfile";
	open my $fh, ">", "$tmpfilename";
	#      print "{{$fasta}}\n";
	print $fh "$overlap_fasta_string"; # "$fasta";
	close $fh;
	#my $support_string = ($do_support)? '': ' -nosupport ';
	my $FT_cl = "FastTree -wag -gamma -bionj $support_string $tmpfilename  2>> FT.log";
	my $newick_string = `$FT_cl`;
	print STDERR "after fasttree \n";
# print STDERR "number of keys of gg_hash: ", scalar keys %$gg_hashref, "\n";
	my $taxonified_newick = taxonify_newick($newick_string, $gg_hashref);
	print STDERR "after taxonify_newick \n";
#	print STDERR "$taxonified_newick \n";
	# reroot it
	# first construct tree obj:
	my $parser = CXGN::Phylo::Parse_newick->new( $taxonified_newick, 1 );
	my $tree = $parser->parse( CXGN::Phylo::BasicTree->new() );

#	find_cycle($parser);

#find_cycle($overlap_obj);
#exit;
	$tree->impose_branch_length_minimum();
	$tree->show_newick_attribute("species");
	$tree->set_show_standard_species(0);
	$tree->get_root()->recursive_implicit_names();
	$tree->get_root()->recursive_implicit_species();

	# if species is empty, (because no [species=x] in newick), get it from the name...
	#$tree->set_missing_species_from_names();

	$tree->set_show_standard_species(1);
	$tree->set_species_standardizer( CXGN::Phylo::Species_name_map->new() );

### make the tree binary
	my @new_root_point;
	{ # if root has > 2 children, reroot to pt on a root-child branch, so new root has 2 children.
	  my @root_children = $tree->get_root()->get_children();
	  if ( scalar @root_children > 2 ) {
            @new_root_point = ( $root_children[0], 0.9 * $root_children[0]->get_branch_length() );
            $tree->reset_root_to_point_on_branch(@new_root_point);
	  }
	}
	# binarify every non-binary node. At present doesn't attempt to choose in a smart way
	# among the various possible resolutions
	$tree->make_binary( $tree->get_root() ); # urec requires binary tree.
### 

	$tree = reroot( $tree, $reroot_method, $species_tree );
	my $rerooted_tree_newick = $tree->generate_newick();



	
	 $tree->decircularize(); # so can be properly garbage collected!
#find_cycle($parser);
#find_cycle($tree);
	print STDERR "after reroot \n";
	print "$rerooted_tree_newick \n\n";
# exit;
      }else{
	print "\n";
      }

      } else {
	print "$string_to_print\n\n";
      }
      $state = 'idline';
    }else{
      $fasta .= $_;
    }
  }
}



sub reroot {
  my $tree             = shift;
  my $reroot_method    = shift || 'mindl';
my $species_tree = shift;
  my $species_tree_arg = shift;	# file name of species tree newick file. overrides default species tree

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

    if(0){
    my $species_tree_newick;
    if ( defined $species_tree_arg ) {
      my $species_tree_file = CXGN::Phylo::File->new($species_tree_arg);
      find_cycle($species_tree_file);
      $species_tree_newick = $species_tree_file->get_tree_string();
    } else {
      $species_tree_newick =
	'( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

      # "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
    }
    my $species_tree_parser = CXGN::Phylo::Parse_newick->new( $species_tree_newick, 1 );
    my $species_tree = $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

    $species_tree->set_missing_species_from_names(); # get species from name if species undef
    $species_tree->impose_branch_length_minimum();
    $species_tree->collapse_tree();
    $species_tree->get_root()->recursive_implicit_names();
    $species_tree->get_root()->recursive_implicit_species();
  }


    my $spec_bit_hash = $tree->get_species_bithash($species_tree);

    $tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);
    $species_tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);

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

sub get_species_tree{
  my $species_tree_arg = shift;	# file name of species tree newick file. default species tree as specified below
  my $species_tree_newick;
  if ( defined $species_tree_arg ) {
    my $species_tree_file = CXGN::Phylo::File->new($species_tree_arg);
    find_cycle($species_tree_file);
    $species_tree_newick = $species_tree_file->get_tree_string();
  } else {
    $species_tree_newick =
      '( chlamydomonas[species=Chlamydomonas_reinhardtii]:1, ( physcomitrella[species=Physcomitrella_patens]:1, ( selaginella[species=Selaginella_moellendorffii]:1, ( loblolly_pine[species=Pinus_taeda]:1, ( amborella[species=Amborella_trichopoda]:1, ( ( date_palm[species=Phoenix_dactylifera]:1, ( ( foxtail_millet[species=Setaria_italica]:1, ( sorghum[species=Sorghum_bicolor]:1, maize[species=Zea_mays]:1 ):1 ):1, ( rice[species=Oryza_sativa]:1, ( brachypodium[species=Brachypodium_distachyon]:1, ( wheat[species=Triticum_aestivum]:1, barley[species=Hordeum_vulgare]:1 ):1 ):1 ):1 ):1 ):1, ( columbine[species=Aquilegia_coerulea]:1, ( ( ( ( ( ( ( ( ( ( tomato[species=Solanum_lycopersicum]:1, potato[species=Solanum_tuberosum]:1 ):1, eggplant[species=Solanum_melongena]:1 ):1, pepper[species=Capsicum_annuum]:1 ):1, tobacco[species=Nicotiana_tabacum]:1 ):1, petunia[species=Petunia]:1 ):1, sweet_potato[species=Ipomoea_batatas]:1 ):1, ( arabica_coffee[species=Coffea_arabica]:1, robusta_coffee[species=Coffea_canephora]:1 ):1 ):1, snapdragon[species=Antirrhinum]:1 ):1, ( ( sunflower[species=Helianthus_annuus]:1, lettuce[species=Lactuca_sativa]:1 ):1, carrot[species=Daucus_carota]:1 ):1 ):1, ( grape[species=Vitis_vinifera]:1, ( ( eucalyptus[species=Eucalyptus_grandis]:1, ( ( orange[species=Citrus_sinensis]:1, clementine[species=Citrus_clementina]:1 ):1, ( ( cacao[species=Theobroma_cacao]:1, cotton[species=Gossypium_raimondii]:1 ):1, ( papaya[species=Carica_papaya]:1, ( (turnip[species=Brassica_rapa]:1, (salt_cress[species=Thellungiella_parvula]:1, (Thellungiella_h[species=Thellungiella_halophila]:0.01, Thellungiella_s[species=Thellungiella_salsuginea]:0.01):1 ):1):1,(red_shepherds_purse[species=Capsella_rubella]:1, ( arabidopsis_thaliana[species=Arabidopsis_thaliana]:1, arabidopsis_lyrata[species=Arabidopsis_lyrata]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1, ( ( ( peanut[species=Arachis_hypogaea]:1, ( ( soy[species=Glycine_max]:1, pigeon_pea[species=Cajanus_cajan]:1 ):1, ( medicago[species=Medicago_truncatula]:1, lotus[species=Lotus_japonicus]:1 ):1 ):1 ):1, ( hemp[species=Cannabis_sativa]:1, ( ( ( apple[species=Malus_domestica]:1, peach[species=Prunus_persica]:1 ):1, woodland_strawberry[species=Fragaria_vesca]:1 ):1, cucumber[species=Cucumis_sativus]:1 ):1 ):1 ):1, ( ( castorbean[species=Ricinus_communis]:1, cassava[species=Manihot_esculenta]:1 ):1, ( poplar[species=Populus_trichocarpa]:1, flax[species=Linum_usitatissimum]:1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 ):1 )';

    # "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
  }
  my $species_tree_parser = CXGN::Phylo::Parse_newick->new( $species_tree_newick, 1 );
  my $species_tree = $species_tree_parser->parse( CXGN::Phylo::BasicTree->new() );

  $species_tree->set_missing_species_from_names(); # get species from name if species undef
  $species_tree->impose_branch_length_minimum();
  $species_tree->collapse_tree();
  $species_tree->get_root()->recursive_implicit_names();
  $species_tree->get_root()->recursive_implicit_species();

  return $species_tree;
}



sub taxonify_newick{
  my $newick = shift;
  my $seqid_species = shift;
  # my %seqid_species = ();
  # if (defined $gg_filename and -f $gg_filename) {
  #   open my $fh_gg, "<", "$gg_filename";
  #   while (<$fh_gg>) {
  #     my @cols = split(" ", $_);
  #     my $species = shift @cols;
  #     $species =~ s/:$//;	# remove final colon if present.
  #     for (@cols) {
  # 	$seqid_species{$_} = $species;
  #     }
  #   }
  # }	      
  # done storing gg_file info in hash %seqid_species
  $newick =~ s/\s+$//;		# remove final whitespace
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
      if (exists $seqid_species->{$seq_id}) {
	my $species = $seqid_species->{$seq_id};
	$seq_id .= '[species=' . $species . ']';
      } else {
# print STDERR "number of keys of gg_hash: ", scalar keys %$seqid_species, "\n";
	warn "sequence $seq_id; no corresponding species found.\n";
      }
      s/ ([\(,]) \s* ([^\(,\):]+) \s* $/$1$seq_id/x;  
    }
  }
  $new_newick = join('', @leaf_lines) . $last_bit;
  return $new_newick;
}

sub store_gg_info{	       # read in gene-genome association file 
  # store in hash (keys: seqids, values: species). e.g. ('ATxxxx' => Arabidopsis_thaliana');
  my $gg_filename = shift;
  my %seqid_species = ();
  if (defined $gg_filename) {
    if ( -f $gg_filename) {
      open my $fh_gg, "<", "$gg_filename";
      while (<$fh_gg>) {
	my @cols = split(" ", $_);
	my $species = shift @cols;
	$species =~ s/:$//;	# remove final colon if present.
	for (@cols) {
	  if (exists $seqid_species{$_}) {
	    warn "key $_ already stored with species: ", $seqid_species{$_}, "\n";
	  } else {
	    $seqid_species{$_} = $species;
	  }
	}
      }
      close $fh_gg;
    }else {     # done storing gg_file info in hash %seqid_species
    die "$gg_filename: no such file.\n";
  }
} else {
  die "gg filename undefined. \n";
}
return \%seqid_species; 
}
