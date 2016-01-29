#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
# use IPC::Run3;
use List::Util qw (min max sum);

# takes as input the output from clumper.pl
#
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
# use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
use lib $libdir;

use Phyml;
use CXGN::Phylo::Overlap;
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';

   my $AM35 = {
               'Amborella_trichopoda' => 1,

               'Phoenix_dactylifera'     => 1, # date palm
                         'Musa_acuminata' => 1,          # banana
                         'Elaeis_guineensis' => 1,       # oil palm

                         'Zea_mays'                => 1, # maize
                         'Panicum_virgatum' => 1,        # switchgrass
                         'Panicum_hallii' => 1, 
                         'Phyllostachys_heterocycla' => 1, # bamboo, AM ??
                         'Brachypodium_distachyon' => 1,
                         'Sorghum_bicolor'         => 1,
                         'Oryza_sativa'            => 1, # rice
                         'Hordeum_vulgare'         => 1, # barley
                         'Triticum_urartu'       => 1, # diploid wheat relative

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
                      }, 


my %species_abcfile = (
      Medicago_truncatula => '/home/tomfy/Aug2015multiquery/medicago_query/sapu2780set/fams200_3/2780set.abc',
      Lotus_japonicus => '/home/tomfy/Aug2015multiquery/lotus_japonicus_query/sapu1940set/fams200_3/1940set.abc',
      Carica_papaya => '/home/tomfy/Aug2015multiquery/papaya_query/sapu19x_1291set/fams200_3/1291set.abc',
      Solanum_lycopersicum => '/home/tomfy/Aug2015multiquery/tomato_query/sapu1133set/fams200_3/1133set.abc',
      Oryza_sativa => '/home/tomfy/Aug2015multiquery/rice_query/sapu1546set/fams200_3/1546set.abc',
      Phoenix_dactylifera => '/home/tomfy/Aug2015multiquery/date_palm_query/sapu879set/fams200_3/879set.abc',
      );

my $clump_qid_filename = undef;
my $gg_filename = undef;
my $abcfilelist = undef;
GetOptions(
      'qids=s' => \$clump_qid_filename,
      'abcfilelist=s' => \$abcfilelist,
	   'gg_file=s'           => \$gg_filename, #
      );
if(defined $abcfilelist){
   open my $fhlist, "<", "$abcfilelist" or warn "couldn't open $abcfilelist for reading. exiting.\n";
   %species_abcfile = ();
   for(<$fhlist>){
      if(/^\s*(\S+)\s+(\S+)\s*$/){
         species_abcfile{$1} = $2;
      }
   }
}
my $gg_hashref = store_gg_info($gg_filename);

   open my $fhin, "<", "$clump_qid_filename" or die "couldn't open $clump_qid_filename for reading.\n";
   my %clumpidnumber_qidset = ();
my %clumpidnumber_clamidset = ();
   my %qid_clumpidnumber = ();
   my $clump_id_number = 1;
   while (my $idline = <$fhin>) {
      next if($idline =~ /^\s*#/);
      my @cols = split(" ", $idline);
      my ($n_qids_clumped, $clump_qids_string, $clump_amcladeids_string) = @cols[3,9,11];
      my @clump_qids = split(",", $clump_qids_string);
      $clump_amcladeids_string =~ s/([|])//g;
      my @clump_amcladeids = split(",", $clump_amcladeids_string);
      $clumpidnumber_qidset{$clump_id_number} = {};
      for my $qid (@clump_qids) {
         $clumpidnumber_qidset{$clump_id_number}->{$qid} = 1;
         if (exists $qid_clumpidnumber{$qid}) {
            warn "query id $qid in > 1 clump? \n";
         } else {
            $qid_clumpidnumber{$qid} = $clump_id_number;
         }

      }                            # end loop over query ids in clump
 $clumpidnumber_clamidset{$clump_id_number} = {};
      for my $clamid (@clump_amcladeids) {
         $clumpidnumber_clamidset{$clump_id_number}->{$clamid} = 1;
      }  
#print STDERR "clump number:  $clump_id_number   qids in clump: ", scalar keys %{$clumpidnumber_qidset{$clump_id_number}}, "\n";
      $clump_id_number++;
   }
   close $fhin;

   my %clumpidnumber_allidset = (); 
   while (my($sp,$abc_filename) = each %species_abcfile) {
      open my $fhabc, "<", "$abc_filename" or die "couldn't open $abc_filename for reading. \n";
      while (my $abcline = <$fhabc>) {
         my @abccols = split(" ", $abcline);
         my ($id1, $id2, $ev) = @abccols[0,1,2];
         if (exists $qid_clumpidnumber{$id1}) { # this qid is part of one of the clumps
            my $clump_id_number = $qid_clumpidnumber{$id1};
            if (exists $clumpidnumber_allidset{$clump_id_number}) {
               $clumpidnumber_allidset{$clump_id_number}->{$id2}++;
            } else {
               $clumpidnumber_allidset{$clump_id_number} = {$id2 => 1};
            }
         } else {
# this query id is not of interest, do nothing with this line
         }
      }
      close $fhabc;
   }

   my @sorted_clump_numbers = sort { scalar keys %{$clumpidnumber_allidset{$b}} <=> scalar keys %{$clumpidnumber_allidset{$a}} } keys %clumpidnumber_allidset;
   print STDERR "# clump id, n_qids, n_all_ids: \n";
   for my $clump_id (@sorted_clump_numbers) {
      print STDERR "$clump_id  ", 
            scalar keys %{$clumpidnumber_qidset{$clump_id}}, "  ", 
            scalar keys %{$clumpidnumber_allidset{$clump_id}}, "\n";

      my @ids =  sort { $clumpidnumber_allidset{$clump_id}->{$b} <=> $clumpidnumber_allidset{$clump_id}->{$a} } keys %{$clumpidnumber_allidset{$clump_id}};
      for (@ids) {
         my $the_sp = $gg_hashref->{$_};
         my $AM = (exists $AM35->{$the_sp})? '1' : '0';
         my $inclade = (exists $clumpidnumber_clamidset{$clump_id}->{$_})? '1' : '0';
         print "clump_$clump_id  $_  ", $clumpidnumber_allidset{$clump_id}->{$_}, " $AM  $inclade  \n";
      }
   }
