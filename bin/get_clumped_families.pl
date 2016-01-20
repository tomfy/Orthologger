#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
# use IPC::Run3;
use List::Util qw (min max sum);

# no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

# use File::Basename 'dirname';
# use Cwd 'abs_path';
# my ( $bindir, $libdir );

# BEGIN {	    # this has to go in Begin block so happens at compile time
#    $bindir =
#      dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
#    $libdir = $bindir . '/../lib';
#    $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
# }
# use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
# use lib $libdir;

# use Phyml;
# use CXGN::Phylo::Overlap;
# use CXGN::Phylo::Parser;
# use CXGN::Phylo::BasicTree;
# use CXGN::Phylo::File;
# use CXGN::Phylo::Species_name_map;

#use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';


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

GetOptions(
           'qids=s' => \$clump_qid_filename,
           #	   'gg_file=s'           => \$gg_filename, #
          );

#my $gg_hashref = store_gg_info($gg_filename);

open my $fhin, "<", "$clump_qid_filename" or die "couldn't open $clump_qid_filename for reading.\n";
my %clumpidnumber_qidset = ();
my %qid_clumpidnumber = ();
my $clump_id_number = 1;
while (my $idline = <$fhin>) {
   my @cols = split(" ", $idline);
   my ($n_qids_clumped, $clump_qids_string) = @cols[3,9];
   my @clump_qids = split(",", $clump_qids_string);
   $clumpidnumber_qidset{$clump_id_number} = {};
   for my $qid (@clump_qids) {
      $clumpidnumber_qidset{$clump_id_number}->{$qid} = 1;
      if (exists $qid_clumpidnumber{$qid}) {
         warn "query id $qid in > 1 clump? \n";
      } else {
         $qid_clumpidnumber{$qid} = $clump_id_number;
      }

   }                            # end loop over query ids in clump
   print STDERR "clump number:  $clump_id_number   qids in clump: ", scalar keys %{$clumpidnumber_qidset{$clump_id_number}}, "\n";
   $clump_id_number++;
}

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

for my $clump_id (@sorted_clump_numbers) {
   print STDERR "$clump_id  ", scalar keys %{$clumpidnumber_allidset{$clump_id}}, "\n";
   my @ids =  sort { $clumpidnumber_allidset{$clump_id}->{$b} <=> $clumpidnumber_allidset{$clump_id}->{$a} } keys %{$clumpidnumber_allidset{$clump_id}};
   for (@ids) {
      print "clump_$clump_id  $_  ", $clumpidnumber_allidset{$clump_id}->{$_}, "\n";
   }
}
