#!/usr/bin/perl -w
use strict;

no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {
    $bindir =
      dirname( abs_path(__FILE__) );    # this has to go in Begin block so happens at compile time
    $libdir = $bindir . '/../lib';
    $libdir = abs_path($libdir);        # collapses the bin/../lib to just lib
}
use lib $libdir;

use CXGN::Phylo::Overlap;

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
my $n_taxa = shift || 21;
my $min_taxa = 4;
my $nongap_fraction = 0.8;
my $state = 'idline';		# other possible values: newick
my ($qid, $fam_size, $taxa, $idline, $newick, $do);
my $support_string = '-nosupport'; # set to '' to get branch support calculation.
while (<>) {
  if ($state eq 'idline') {
    my @cols = split(" ", $_);
    ($qid, $fam_size, $taxa, $do) = @cols[1,4,5,6];
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
    $state = 'newick';
    $newick = '';
  } elsif ($state eq 'newick') {
  
    if (/^\s*$/) { # blank line after sequence -> process the sequence.
  $newick .= $_;
      chomp $idline;
  my @fields = split(" ", $idline);
  my $id1 = $fields[1]; # Id Medtr...  family
      print "$idline  $do \n";
      if ($do) {



	open my $fh, ">", "ZZZZZtmp";
	#      print "{{$fasta}}\n";
	print $fh "$newick";
	close $fh;
	#my $support_string = ($do_support)? '': ' -nosupport ';
	my $rp_cl = "reroot_prune.pl --input ZZZZZtmp  --sequence $id1 --ggfile $gg_filename --prune_threshold 3  ";
	my $newick_string = `$FT_cl`;
	print "$newick_string \n";
      } else {
	print "\n";
      }
      $state = 'idline';
    }
  }
}
