#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);
use Graph;

# clustering based on (approximate) reciprocal best matches
# input is abc file (i.e. each line has id1 id2 e-value)
# make a graph with edges joining reciprocal best matches
# (or approximate rbms - e-value in both directions is withing $F of best)
# cluster are the biconnected components of the graph.

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;
use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';

use constant LOG10 => log(10.0);
use constant NOTANID => '_not_an_actual_id_';
my $logF = 1.0;
my $max_matches = 600;
my $max_fam_size = 400;
my $min_ev = 1e-180;
my $max_sim = 181.0;
my $min_sim = 0;
my $weed = 1;
my $clusters_filename = 'clusters';
my $singles_filename = 'singles';
# store best match of each species

my $abc_filename = shift;
my $gg_filename = shift || '/home/tomfy/Aug2015multispeciesquery/55set.gg';
my $geneid_sp = store_gg_info("$gg_filename");
print STDERR "\n", "# Done storing gg info.\n";
my %species_count =  ();
my %id1_id2ev = ();
my %id_selfev = ();
my $id2_ev = {};
my %singles = ();
my @clusters = ();
my %clusteridn_ids = ();
my %id_clusters = ();
my %id_ref = (); # key is an id (string), value is ref to scalar containing that string
my %id_number = ();

my $cluster_id_number = 1;
while ( my $line = <>) {
   $line =~ /^\s*(\S+)/;
   my @ids = split(",", $1);
   if (scalar @ids == 1) {
      $singles{$ids[0]} = 1;
   } elsif (scalar @ids > 1) {
      push @clusters, \@ids;
      %clusteridn_ids{$cluster_id_number} = \@ids;
      for (@ids) {
         if (exists $id_clusters{$_}) {
            push @{$id_clusters{$_}}, $cluster_id_number;
         } else {
            $id_clusters{$_} = [$cluster_id_number];
         }
      }
      $cluster_id_number++;
   } else {
      warn "Input line should be , separated list of ids, is: $line \n";
      next;
   }
}
### done storing query ids of clusters.

### read in blast match info from abc file:
my $id_counter = 0;
open my $fh_abc, "<", "$abc_filename" or die "Couldn't open $abc_filename for reading. Exiting. \n";
my $old_id1 = NOTANID;
my $count_matches = 0;
while (my $line = <$fh_abc>) {
   my ($id1, $id2, $ev) = split(" ", $line);
   if ($id1 ne $old_id1 and ($old_id1 ne NOTANID) ) {
      # if $old_id1 is a single, output its family:
      if (exists $singles{$old_id1}) {
         print "$old_id1 \n";
         my @sorted_id2s = sort {$id2_ev->{$a} <=> $id2_ev->{$b}} keys %$id2_ev;
       #  while (my ($id2, $ev) = each %$id2_ev) {
            for my $id2 (@sorted_id2s){
               my $ev = $id2_ev->{$id2};
            printf ("  %s %7.2f \n", $id2, ev_to_sim($ev));
         }
      } else {
         $id1_id2ev{$old_id1} = $id2_ev;
      }
      $count_matches = 0;
      #   $id1_id2ev{$id1} = {};
      $id2_ev = {};
      %species_count = ();
   }
   if(!exists $id_number{$id1}){
      $id_number{$id1} = $id_counter;
      $id_counter++;
   }
 if(!exists $id_number{$id2}){
      $id_number{$id2} = $id_counter;
      $id_counter++;
   }

   $count_matches++;
   next if($count_matches >= $max_matches);
   #   $id1_id2ev{$id1}->{$id2} = $ev;
   $id2_ev->{$id2} = $ev;
   my $self_ev = undef;
   if ($id1 eq $id2) {
      $self_ev = $ev;
      $id_selfev{$id1} = $self_ev;
      next;
   }
   my $sp2 = $geneid_sp->{$id2};
   $species_count{$sp2}++; # skip if not one of the specified species.
   $old_id1 = $id1;
}
print STDERR "# Done reading in abc data.\n";
####### Done reading in abc data #####


for my $cluster (@clusters) {
   my @cluster_qids = @$cluster;
   my %id2_evs = ();
   my $cl_size = scalar @ids;
   my $max_ev = -1;
   for my $id1 (@cluster_qids) {
      while (my ($id2, $ev) = each %{$id1_id2ev{$id1}}) {
         if (!exists $id2_evs{$id2}) {
            $id2_evs{$id2} = [$ev];
         } else {
            push @{$id2_evs{$id2}}, $ev;
         }
         $max_ev = $ev if($ev > $max_ev);
      }
      $min_sim = ($max_ev == 0)? $max_sim : -log($max_ev)/LOG10;
   }
   my %id2_avgsim = ();
   while (my ($id2, $evals) = each %id2_evs) {
      my $sim_sum = 0;
      for (@$evals) {
         $sim_sum += ($_ == 0)? $max_sim : -log($_)/LOG10;
      }
      my $n_evs = scalar @$evals;
      $sim_sum += ($cl_size - $n_evs)*$min_sim;
      my $avg_sim = $sim_sum/$cl_size;
      $id2_avgsim{$id2} = $avg_sim;
   }
   my @sid2s = sort { $id2_avgsim{$b} <=> $id2_avgsim{$a} } keys %id2_avgsim;
   print "$cluster_ids \n";
   my $fam_size = scalar @sid2s;
   $fam_size = $max_fam_size if($fam_size > $max_fam_size);
   for (0..$fam_size-1) {
      my $id2 = $sid2s[$_];
      printf ("  %s %7.2f \n", $id2, $id2_avgsim{$id2});
   }
}

####### 
while ( my $line = <>) {
   chomp $line;
   my %id2_evs = ();
   $line =~ /^\s*(\S+)/;
   my $cluster_ids = $1;
   my @ids = split(",", $cluster_ids);
   my $cl_size = scalar @ids;
   my $max_ev = -1;
   for my $id1 (@ids) {
      while (my ($id2, $ev) = each %{$id1_id2ev{$id1}}) {
         if (!exists $id2_evs{$id2}) {
            $id2_evs{$id2} = [$ev];
         } else {
            push @{$id2_evs{$id2}}, $ev;
         }
         $max_ev = $ev if($ev > $max_ev);
      }
      $min_sim = ($max_ev == 0)? $max_sim : -log($max_ev)/LOG10;
   }
   my %id2_avgsim = ();
   while (my ($id2, $evals) = each %id2_evs) {
      my $sim_sum = 0;
      for (@$evals) {
         $sim_sum += ($_ == 0)? $max_sim : -log($_)/LOG10;
      }
      my $n_evs = scalar @$evals;
      $sim_sum += ($cl_size - $n_evs)*$min_sim;
      my $avg_sim = $sim_sum/$cl_size;
      $id2_avgsim{$id2} = $avg_sim;
   }
   my @sid2s = sort { $id2_avgsim{$b} <=> $id2_avgsim{$a} } keys %id2_avgsim;
   print "$cluster_ids \n";
   my $fam_size = scalar @sid2s;
   $fam_size = $max_fam_size if($fam_size > $max_fam_size);
      for(0..$fam_size-1){
         my $id2 = $sid2s[$_];
      printf ("  %s %7.2f \n", $id2, $id2_avgsim{$id2});
   }
}

sub ev_to_sim{
   my $evalue = shift;
   my $max_sim = 181.0;
   my $sim = ($evalue == 0)? $max_sim : -log($evalue)/LOG10;
   return $sim;
}
