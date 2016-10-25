#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);
# use Devel::Size qw (total_size size);
# use Time::HiRes qw (usleep);
use Graph;

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {     # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir); # collapses the bin/../lib to just lib
}
use lib $libdir;

use ClMatches;

# clustering based on (approximate) reciprocal best matches
# input is abc file (i.e. each line has id1 id2 e-value)
# make a graph with edges joining reciprocal best matches
# (or approximate rbms - e-value in both directions is withing $F of best)
# cluster are the biconnected components of the graph.

use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';

use constant LOG10 => log(10.0);
use constant NOTANID => '_not_an_actual_id_';
my $clusters_printed = 0;
my $max_fam_size = 400;
my $max_sim = 181.0;
my $min_sim = 0;
my $clusters_filename = 'clusters';
# store best match of each species

my $abc_filename = shift;

my $gg_filename = shift || '/home/tomfy/Aug2015multispeciesquery/55set.gg';
# my $geneid_sp = store_gg_info("$gg_filename");
# print STDERR "\n", "# Done storing gg info.\n";

my %qid_clusterobjs = ();
my $count_cluster_objs = 0;
my $cluster_id_number = 1;
while ( my $line = <>) {
   if ($line =~ /^\s*(\S+)/) {
      my $clm_obj = ClMatches->new($1, $min_sim, $max_sim); $count_cluster_objs++;
      my @qids = split(",", $1);
      for (@qids) {
         if (!exists $qid_clusterobjs{$_}) {
            $qid_clusterobjs{$_} = [];
         }
         push @{$qid_clusterobjs{$_}}, $clm_obj;
         $cluster_id_number++;
      }
   } else {
      warn "Input line should be , separated list of ids, is: $line \n";
      next;
   }
}

print STDERR "number of cluster objects: ", $count_cluster_objs, "    ", scalar keys %qid_clusterobjs, "\n";
#exit;


### read in blast match info from abc file:
my $id_counter = 0;
open my $fh_abc, "<", "$abc_filename" or die "Couldn't open $abc_filename for reading. Exiting. \n";
my $old_id1 = NOTANID;
my $count_matches = 0;
my $matches_string = '';
my %id2_ev = ();
my $count_lines_read = 0;
my $count_qids_processed = 0;
while (my $line = <$fh_abc>) {
   $count_lines_read++; print STDERR "abc lines read: $count_lines_read \n" if($count_lines_read % 3000000 == 0);
   my ($id1, $id2, $ev) = split(" ", $line);
   if (($id1 ne $old_id1) and ($old_id1 ne NOTANID)) {
      if (exists $qid_clusterobjs{$old_id1}) { # this qid belongs to one of the clusters, so process its matches:
         process_matches(\%qid_clusterobjs, $old_id1, $matches_string, $max_fam_size);
         $count_qids_processed++;
      }
      $matches_string = $line;
      %id2_ev = ();
   } else {
      next if(exists $id2_ev{$id2});
      $id2_ev{$id2} = 1;
      $matches_string .= $line;
   }
   $old_id1 = $id1;
}
if (exists $qid_clusterobjs{$old_id1}) {
   process_matches(\%qid_clusterobjs, $old_id1, $matches_string, $max_fam_size);
   $count_qids_processed++;
}
if($count_qids_processed % 1000000 == 0){
  my %clustobjs = ();
         while(my($k, $v) = each %qid_clusterobjs){
            for(@$v){
               $clustobjs{$_->{clusterqids_str}}++;
            }
         }
         print STDERR "Xnumber of cluster objects remaining: ", scalar keys %clustobjs, "\n";
}
$matches_string = '';
print STDERR "# Done reading in abc data.\n";
print STDERR "# qids in clusters processed: ", $count_qids_processed, "\n";
print STDERR "# qids in clusters yet to do: ", scalar keys %qid_clusterobjs, "\n";
while(my($k, $v) = each %qid_clusterobjs){
   print STDERR $k, "    ", scalar @$v,  "\n";
}
####### Done reading in abc data #####


sub process_matches{
   my $qid_clobjs = shift;
   my $qid = shift;
   my $matches_str = shift;
   my $max_famsize = shift;
   my @the_objs = @{$qid_clobjs->{$qid}};
   for my $the_obj (@the_objs) {
      $the_obj->add_matches($matches_str);
      if ($the_obj->all_qids_done()) {
         print $the_obj->get_qids_str(), "\n", $the_obj->top_n_by_avg_sim($max_famsize);
         $clusters_printed++;
      } else {
      }
   }
   delete $qid_clobjs->{$qid};
}
