#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
use IPC::Run3;

no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {
  $bindir =
    dirname( abs_path(__FILE__) )
      ;	    # this has to go in Begin block so happens at compile time
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
my $gg_filename; # required  cl parameter 
my $do_nj = 0;
my $do_ml = 1;
my $n_bs = 0; # number of bs replicates
my $ml_bs = 0; # by default, do not do ML for bs
# print "[$gg_filename] [$do_nj] [$do_ml] [$n_bs] [$ml_bs] \n";
# my $n_taxa                   = shift || 21;
my $species_tree_newick_file = undef;
my $reroot_method            = 'mindl';
my $nongap_fraction = 0.8;
my $min_overlap_length = 40;
my $nosupport = ''; # default is to do support
GetOptions(
    'gg_file=s'           => \$gg_filename,
    'nj!'          => \$do_nj,
    'ml!'      => \$do_ml,
    'nosupport'          => \$nosupport,
    'n_bs=i' => \$n_bs,
    'ml_bs!' => \$ml_bs, # boolean to do ML bs or not (default: 0)
	   'species_tree=s' => \$species_tree_newick_file,
	   'reroot_method=s' => \$reroot_method, 
	   'nongap_fraction=i' => \$nongap_fraction,
	   'min_overlap_length=i' => \$min_overlap_length
);
if(!defined $gg_filename  or  ! -f $gg_filename){
  die "No gene-genome association file specified. Exiting. \n";
}
print STDERR "gg_file: ", (defined $gg_filename)? $gg_filename: 'undefined', " \n", 
  "do_nj: $do_nj \n",
  "do_ml: $do_ml \n", 
  "branch support: ", $nosupport? '0' : '1', " \n",
  "n_bs: $n_bs \n", "ml_bs: $ml_bs \n", 
  "species_tree: ", (defined $species_tree_newick_file)? $species_tree_newick_file : 'undefined', " \n",
  "reroot_method: $reroot_method \n",
  "nongap_fraction: $nongap_fraction \n",
  "min_overlap_length: $min_overlap_length \n";
#my $min_overlap_length = 40;
# my $min_taxa           = 4;
#my $nongap_fraction    = 0.8;
my $state              = 'idline'; # other possible values: fasta
my ( $qid, $fam_size, $taxa, $idline, $fasta);
my $support_string = ($nosupport)? '-nosupport' : '';  # controls whether FastTree calculates local branch support numbers.
my $gg_hashref = store_gg_info($gg_filename);

#print STDERR "number of keys of gg_hash: ", scalar keys %$gg_hashref, "\n";

my $species_tree = get_species_tree($species_tree_newick_file)
  ; # get species tree from file (or use default if file undef or doesn't exist);
my $n_ml = 0;
while (<>) {
#  print STDERR "state: $state; line read: $_";
  if ( $state eq 'idline' ) {
    my @cols = split( " ", $_ );
    ( $qid, $fam_size, $taxa) = @cols[ 1, 4, 5];
 #   print STDERR "$qid, $fam_size\n";

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
    $state  = 'fasta';
    $fasta  = '';
  } elsif ( $state eq 'fasta' ) {

    #    $fasta .= $_;
    if (/^\s*$/) { # blank line after sequence -> process the sequence.
      chomp $idline;
      my $string_to_print = "$idline   ";
      my $rng_seed = srand();
      if ($fasta ne '') {

	#	print STDERR "[[[$fasta]]]\n";
#	print STDERR "before overlap \n";
	my $overlap_obj =
	  CXGN::Phylo::Overlap->new( $fasta, $nongap_fraction,
				     $rng_seed );
	my $overlap_fasta_string =
	  $overlap_obj->overlap_fasta_string('');
	my $overlap_length = $overlap_obj->get_overlap_length();
#	print STDERR  "after overlap. overlap length: $overlap_length.\n";

#	print STDERR "[[[$overlap_fasta_string]]]\n";
	$string_to_print .= "$overlap_length ";
#	print "$string_to_print \n";
	if ( $overlap_length >= $min_overlap_length ) {
	  ################## do actual data - NJ  ########################################
	  if($do_nj){
#print STDERR $overlap_fasta_string;
	  my $nj_newick = run_quicktree($overlap_fasta_string);
#print STDERR "after run_quicktree \n"; exit;
#print STDERR "$nj_newick \n";
	  my $taxonified_nj_newick =
	    taxonify_newick( $nj_newick, $gg_hashref );
#print STDERR "after taxonify_newick \n";
	  my $rr_nj_newick =
	    newick2rerootednewick( $taxonified_nj_newick,
				   $reroot_method, $species_tree );
#print STDERR "after newick2rerootednewick \n";
	print "$string_to_print \n";
	  print "NJ  $rr_nj_newick \n\n";
#exit;
	}
	  ################## do actual data - ML  ########################################
	  if($do_ml){
	#    print "$do_ml exiting\n"; exit;
#	  print STDERR "do_ml: [$do_ml]. \n", $overlap_fasta_string, "\n";
print STDERR "fasttree cl: FastTree -wag -gamma -bionj $support_string \n";
	  my ($ml_newick, $ml_stderr) = run_fasttree($overlap_fasta_string, "FastTree -wag -gamma -bionj $support_string ");
# exit;
#  'FastTree -wag -gamma -bionj -nosupport
# print STDERR "ML newick: ", $ml_newick, "\n"; #exit;
	  my $taxonified_ml_newick =
	    taxonify_newick( $ml_newick, $gg_hashref );
             #       print STDERR  "after taxonify_newick \n";
#print "taxonified newick: $taxonified_ml_newick \n XXX \n";
	  #    print STDERR "before ml newick2rerootednewick. taxonified newick: \n", $taxonified_ml_newick, "\n";

	  my $rr_ml_newick =
	    newick2rerootednewick( $taxonified_ml_newick,
				   $reroot_method, $species_tree );
#	      print STDERR "after ml newick2rerootednewick \n";
# exit;
$n_ml++;
	print "$string_to_print \n";
	  print "ML  $rr_ml_newick \n\n"; # exit;
}
#exit;
	  ######################## do bootstraps #########################################
	  ###### NJ ########
	  for ( my $i = 1 ; $i <= $n_bs ; $i++ ) {
	    my $bs_overlap_fasta_string =
	      $overlap_obj->bootstrap_overlap_fasta_string('');
	    my $nj_newick = run_quicktree($bs_overlap_fasta_string);
	    my $taxonified_nj_newick =
	      taxonify_newick( $nj_newick, $gg_hashref );
	#      print STDERR "before nj newick2rerootednewick $taxonified_nj_newick   XXX\n";

	    my $rr_nj_newick =
	      newick2rerootednewick( $taxonified_nj_newick,
				     $reroot_method, $species_tree );
#	      print STDERR "after nj newick2rerootednewick \n";

	print "$string_to_print \n";
	    print "NJ_BS  $rr_nj_newick \n\n";
	    ######### ML ########
	    if ($do_ml and $ml_bs) {
	      my ($ml_newick, $ml_stderr) = run_fasttree($bs_overlap_fasta_string, "FastTree -wag -gamma -bionj $support_string");
	      my $taxonified_ml_newick =
		taxonify_newick( $ml_newick, $gg_hashref );
	      my $rr_ml_newick =
		newick2rerootednewick( $taxonified_ml_newick,
				       $reroot_method, $species_tree );
#	      print STDERR "after ml newick2rerootednewick \n";

	print "$string_to_print \n";
	      print "ML_BS  $rr_ml_newick \n\n";
	    }
	  }			# loop over bootstraps
	}			# overlap >= minlength
	else {
	  print "\n";
	}	# this query has no alignment, just print a blank line
      }		# if ($do) i.e. there is an alignment for this query
      else {	# 
	print "$string_to_print\n\n";
      }
      $state = 'idline';
    } else {			# not an empty line - append to $fasta
      $fasta .= $_;
    }
  }
}


sub get_species_tree {
  my $species_tree_arg = shift
    ; # file name of species tree newick file. default species tree as specified below
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

  return $species_tree;
}

sub taxonify_newick {
  my $newick        = shift;
  my $seqid_species = shift;

# print "newick: \n [$newick] \n";
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
  $new_newick =~ s/ ([\(,]) \s* ([^,\):]+) \s* ([,\):])/$1$2\n$3/xg
    ;				# add newline after each leaf
  my @leaf_lines = split( "\n", $new_newick );
  my $last_bit = pop @leaf_lines;
  for (@leaf_lines) {
    if (/\[species=.*\]/) { # species info already there for this leaf - leave it alone
    } else {

      #      print "XXX: $_ \n";
      / [\(,] \s* ([^\(,\):]+) \s* $/x;
      my $seq_id = $1;

      #      print "YYY: $seq_id \n";
      if ( exists $seqid_species->{$seq_id} ) {
	my $species = $seqid_species->{$seq_id};
	$seq_id .= '[species=' . $species . ']';
      } else {

	# print STDERR "number of keys of gg_hash: ", scalar keys %$seqid_species, "\n";
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
  if ( defined $gg_filename ) {
    if ( -f $gg_filename ) {
      open my $fh_gg, "<", "$gg_filename";
      while (<$fh_gg>) {
	my @cols = split( " ", $_ );
	my $species = shift @cols;
	$species =~ s/:$//;	# remove final colon if present.
	for (@cols) {
	  if ( exists $seqid_species{$_} ) {
	    warn "key $_ already stored with species: ",
	      $seqid_species{$_}, "\n";
	  } else {
	    $seqid_species{$_} = $species;
	  }
	}
      }
      close $fh_gg;
    } else {	    # done storing gg_file info in hash %seqid_species
      die "$gg_filename: no such file.\n";
    }
  } else {
    die "gg filename undefined. \n";
  }
  return \%seqid_species;
}

# sub ovrlp2njtree{
#   my $overlap_string = shift;
#   my $gg_hashref = shift;
# 	  my $tmpfilename = "PID" . $$ . "_tmpfile";
# 	open my $fh, ">", "$tmpfilename";
# 	#      print "{{$fasta}}\n";
# 	print $fh "$overlap_string"; # "$fasta";
# 	close $fh;
# 	my $support_string = (0)? '': ' -nosupport ';
# 	my $FT_cl = "FastTree -wag -gamma -bionj $support_string $tmpfilename  2>> FT.log";
# 	my $newick_string = run_fasttree($overlap_string, $FT_cl );  #`$FT_cl`;
# 	print STDERR "after fasttree \n";
# # print STDERR "number of keys of gg_hash: ", scalar keys %$gg_hashref, "\n";
# 	my $taxonified_newick = taxonify_newick($newick_string, $gg_hashref);
# 	print STDERR "after taxonify_newick \n";
# #	print STDERR "$taxonified_newick \n";
# 	# reroot it
# 	# first construct tree obj:
# 	my $parser = CXGN::Phylo::Parse_newick->new( $taxonified_newick, 1 );
# 	my $tree = $parser->parse( CXGN::Phylo::BasicTree->new() );

# #	find_cycle($parser);

# #find_cycle($overlap_obj);
# #exit;
# 	$tree->impose_branch_length_minimum();
# 	$tree->show_newick_attribute("species");
# 	$tree->set_show_standard_species(0);
# 	$tree->get_root()->recursive_implicit_names();
# 	$tree->get_root()->recursive_implicit_species();

# 	# if species is empty, (because no [species=x] in newick), get it from the name...
# 	#$tree->set_missing_species_from_names();

# 	$tree->set_show_standard_species(1);
# 	$tree->set_species_standardizer( CXGN::Phylo::Species_name_map->new() );

# ### make the tree binary
# 	my @new_root_point;
# 	{ # if root has > 2 children, reroot to pt on a root-child branch, so new root has 2 children.
# 	  my @root_children = $tree->get_root()->get_children();
# 	  if ( scalar @root_children > 2 ) {
#             @new_root_point = ( $root_children[0], 0.9 * $root_children[0]->get_branch_length() );
#             $tree->reset_root_to_point_on_branch(@new_root_point);
# 	  }
# 	}
# 	# binarify every non-binary node. At present doesn't attempt to choose in a smart way
# 	# among the various possible resolutions
# 	$tree->make_binary( $tree->get_root() ); # urec requires binary tree.

# 	$tree = reroot( $tree, $reroot_method, $species_tree );
# 	my $rerooted_tree_newick = $tree->generate_newick();

# 	 $tree->decircularize(); # so can be properly garbage collected!
# #find_cycle($parser);
# #find_cycle($tree);

# 	print STDERR "after reroot \n";
# 	print "$rerooted_tree_newick \n\n";
#   return $rerooted_tree_newick;
# }
# exit;

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

sub run_quicktree {
  my $overlap_fasta_string       = shift;
  my $correction                 = shift || 'kimura';
  my $tmp_overlap_fasta_filename = "PID" . $$ . "_tmp_overlap_fasta";
  open my $fhtmp, ">", "$tmp_overlap_fasta_filename";
  print $fhtmp $overlap_fasta_string, "\n";
  close $fhtmp;

  my $tmp_overlap_stockholm_filename = "PID" . $$ . "_tmp_overlap_stockholm";
  system
    "sreformat stockholm $tmp_overlap_fasta_filename > $tmp_overlap_stockholm_filename";
  my $newick_out;
  if ( $correction eq 'kimura' ) {
    $newick_out = `quicktree -kimura $tmp_overlap_stockholm_filename`;
  } else {
    $newick_out = `quicktree  $tmp_overlap_stockholm_filename`;
  }

  $newick_out =~ s/\s+//g;	# remove whitespace
  $newick_out =~ s/;$//;
  return $newick_out;
}

sub run_fasttree {
  my $overlap_fasta_string = shift;
  my $fasttree_cl          = shift;
  my $intree               = shift;
  if ($intree) {
    open my $fh, ">", "tmp_intree";
    print $fh $intree, "\n";
    close $fh;
    $fasttree_cl .= " -intree tmp_intree ";
  }

#print STDERR "fasttree cl: $fasttree_cl \n"; exit;
  my $fasttree_newick_out = "ft_newick_default_output";
  my $fasttree_stderr_out = "ft_stderr_default_output";
#  print STDERR "fasttree cl: ", $fasttree_cl, "\n";
#  print "AAAAA: ft newickout, stderrout: $fasttree_newick_out, $fasttree_stderr_out \n";
  if (0) { # doesn't seem to work now, did before. Why???
    my $run3_return_value = run3(
				 "$fasttree_cl",        \$overlap_fasta_string,
				 \$fasttree_newick_out, \$fasttree_stderr_out
				);
#    print STDERR "in run_fasttree, run3 return value: [", $run3_return_value, "]\n";
  } else {			# do using a file
    my $temp_file = "PID" . $$ . "_ft_in_tmpfile";
    my $stderr_outfile  = "PID" . $$ . ".stderr";
    open my $fh, ">", "$temp_file";
    #  open my $fhstderr, ">", "$stderr_outfile";
    print $fh $overlap_fasta_string, "\n"; close $fh;
# print STDERR "$fasttree_cl  $temp_file \n"; #exit;
$fasttree_newick_out = `$fasttree_cl $temp_file 2>> $stderr_outfile`;
}


#print STDERR "FT OUT NEWICK: \n[ ", $fasttree_newick_out, " ]\n\n"; exit;
#print STDERR "FT stderr out: \n", $fasttree_newick_out, "\n";
#print STDERR "FT OUT stderr: \n", $fasttree_stderr_out, "\n\n";


  # Gamma(20) LogLk = -5372.474 alpha = 1.562 rescaling lengths by 1.044   # parse ll out of ft stderr output.
  my $fasttree_loglikelihood =
    ( $fasttree_stderr_out =~
      / Gamma [(] \d+ [)] \s+ LogLk \s+ = \s+ ([-] \d+ [.] \d*) \s+ alpha/xm
    ) ? $1 : undef;
  $fasttree_newick_out =~ s/\s+//g;
  $fasttree_newick_out =~ s/;$//;
  return ( $fasttree_newick_out, $fasttree_loglikelihood );
}

sub newick2rerootednewick
  { # take a taxonified newick expression, make tree, reroot is by mindl wrt $species_tree

    # return newick expression
    my $taxonified_newick = shift;
$taxonified_newick = put_min_branchlengths_into_newick($taxonified_newick, 0.01);
# print STDERR "ZZZ: $taxonified_newick \n";
    my $reroot_method     = shift;
    my $species_tree      = shift;
    my $parser = CXGN::Phylo::Parse_newick->new( $taxonified_newick, 0 );
 #   print STDERR "after constructing parser \n";
    my $tree   = $parser->parse( CXGN::Phylo::BasicTree->new() );
    if(!defined $tree){ return 'warning: Couldnt parse taxonified tree in newick2rerootnewick.'; }
    $tree->impose_branch_length_minimum();
    $tree->show_newick_attribute("species");
# print STDERR "AAAA:  ", $tree->generate_newick(), "\n";
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
	@new_root_point = (
			   $root_children[0], 0.9 * $root_children[0]->get_branch_length()
			  );
	$tree->reset_root_to_point_on_branch(@new_root_point);
      }
    }

    # binarify every non-binary node. At present doesn't attempt to choose in a smart way
    # among the various possible resolutions
    $tree->make_binary( $tree->get_root() ); # urec requires binary tree.
#print STDERR "Before reroot.\n";
#print STDERR "newick: \n", $tree->generate_newick(), "\n\n";
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
