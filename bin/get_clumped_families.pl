#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
# use IPC::Run3;
use List::Util qw (min max sum);

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
use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
use lib $libdir;

use Phyml;
use CXGN::Phylo::Overlap;
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';


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
	   'gg_file=s'           => \$gg_filename, #
    );

my $gg_hashref = store_gg_info($gg_filename);

open my $fhin, "<", "$clump_qid_filename" or die "couldn't open $clump_qid_filename for reading.\n";

my $clump_id_number = 0;
while(my $idline = <$fhin>){
    my @cols = split(" ", $idline);
    my ($N, $clump_qids_string) = @cols[3,9];
    my @clump_qids = split(",", $clump_qids_string);
    
    my %clump_ids = ();
    for my $qid (@clump_qids){
	my $qid_species = $gg_hashref->{$qid};
	my $abc_filename = $species_abcfile{$qid_species};
	open my $fhabc, "<", "$abc_filename" or die "couldn't open $abc_filename for reading. \n";
	while(my $abcline = <$fhabc>){
	    my @abccols = split(" ", $abcline); 
	    my ($id1, $id2) = @abccols[0,1];
	 #   print "$qid $id1 $id2 \n";
	    if($id1 eq $qid){
		$clump_ids{$id2}++;
	    }
	}
	close $fhabc;
    } # end loop over query ids in clump
    $clump_id_number++;
    print "clump number:  $clump_id_number   ids in union-clump: ", scalar keys %clump_ids, "\n";

}	  
