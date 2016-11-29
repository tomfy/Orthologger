#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);
use Getopt::Long;
use constant NOTANID => '_not_an_actual_id_';

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

# read clustered qids (output from rbmcl.pl) and create cluster object for each cluster,
# and for each qid store a list of clusters it belongs to (most only belong to 1, but some 
# belong to 2 or more).
# Then for each query id, store matches (id2 and ev) in the cluster obj of each
# cluster the qid is in.
# For each cluster, get the best matches, defined as matches with best average similarity,
# averaged over the qids in cluster

#use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';

my $clusters_printed = 0;
my $max_fam_size = 400;
my $max_sim = 181.0;
my $min_sim = 0;
#my $clusters_filename = 'clusters';
# store best match of each species

my $input_filename = undef;
my $clusters_filename = undef;
my $output_filename = undef;

GetOptions(
	   'input_filename=s'           => \$input_filename,  #
	   'clusters_filename=s'      => \$clusters_filename, # 
           'output_filename=s' => \$output_filename,
           'max_fam_size=i' => \$max_fam_size,
          );

my %qid_clusterobjs = ();
my $count_cluster_objs = 0;
my $cluster_id_number = 1;
open my $fh_clust, "<", "$clusters_filename" or die "couldn't open $clusters_filename for reading. \n";
while ( my $line = <$fh_clust>) {
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
close $fh_clust;

print STDERR "number of cluster objects: ", $count_cluster_objs, "    ", scalar keys %qid_clusterobjs, "\n";

my $count_qids_processed = 0;
### read in blast match info from abc file:
if (file_format_is_abc($input_filename)) {
   open my $fh_in, "<", "$input_filename" or die "Couldn't open $input_filename for reading. Exiting. \n";
   my $old_id1 = NOTANID;
   my $count_matches = 0;
   my $matches_string = '';
   my %id2_ev = ();
   my $count_lines_read = 0;
   $count_qids_processed = 0;
   open my $fh_output, ">", "$output_filename" or die "Couldn't open $output_filename for writing.\n";
   while (my $line = <$fh_in>) {
      $count_lines_read++; print STDERR "abc lines read: $count_lines_read \n" if($count_lines_read % 3000000 == 0);
      my ($id1, $id2, $ev) = split(" ", $line);
      if (($id1 ne $old_id1) and ($old_id1 ne NOTANID)) {
         if (exists $qid_clusterobjs{$old_id1}) { # this qid belongs to one of the clusters, so process its matches:
            process_matches(\%qid_clusterobjs, $old_id1, $matches_string, $max_fam_size, $fh_output);
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
      process_matches(\%qid_clusterobjs, $old_id1, $matches_string, $max_fam_size, $fh_output);
      $count_qids_processed++;
   }
} elsif (file_format_is_iie($input_filename)) { # iie
   open my $fh_in, "<", "$input_filename" or die "Couldn't open $input_filename for reading. Exiting. \n";
   my $old_id1 = NOTANID;
   my $count_matches = 0;
   my $matches_string = '';
   my %id2_ev = ();
   my $count_lines_read = 0;
   $count_qids_processed = 0;
   my $id1;
   open my $fh_output, ">", "$output_filename" or die "Couldn't open $output_filename for writing.\n";
   while (my $line = <$fh_in>) {
      $count_lines_read++; print STDERR "iie lines read: $count_lines_read \n" if($count_lines_read % 3000000 == 0);
      if ($line =~ /^(\S+)/) {  # first line of family
         $id1 = $1;
         #  my ($id1, $id2, $ev) = split(" ", $line);
         if ($old_id1 ne NOTANID) {
            if (exists $qid_clusterobjs{$old_id1}) { # this qid belongs to one of the clusters, so process its matches:
               process_matches(\%qid_clusterobjs, $old_id1, $matches_string, $max_fam_size, $fh_output);
               $count_qids_processed++;
            }
            $matches_string = ''; # $line;
            %id2_ev = ();
         }
      } elsif ($line =~ /^\s+(\S+)\s+(\S+)/) {
         my ($id2, $ev) = ($1, $2);
         next if(exists $id2_ev{$id2});
         $id2_ev{$id2} = 1;
         $matches_string .= "$id1 $line";
      }
      $old_id1 = $id1;
   }
   if (exists $qid_clusterobjs{$old_id1}) {
      process_matches(\%qid_clusterobjs, $old_id1, $matches_string, $max_fam_size, $fh_output);
      $count_qids_processed++;
   }
} else {
   die "input file $input_filename has unrecognized format.\n";
}

print STDERR "# Done reading in blast match data.\n";
print STDERR "# qids in clusters processed: ", $count_qids_processed, "\n";
print STDERR "# qids in clusters yet to do: ", scalar keys %qid_clusterobjs, "\n";

####### Done reading in blast match data (abc or iie) #####


sub process_matches{
   my $qid_clobjs = shift;
   my $qid = shift;
   my $matches_str = shift;
   my $max_famsize = shift;
   my $fh_out = shift;
   my @the_objs = @{$qid_clobjs->{$qid}};
   for my $the_obj (@the_objs) {
      $the_obj->add_matches($matches_str);
      if ($the_obj->all_qids_done()) {
         print $fh_out $the_obj->get_qids_str(), "\n", $the_obj->top_n_by_avg_sim($max_famsize);
         $clusters_printed++;
      } else {
      }
   }
   delete $qid_clobjs->{$qid};
}

sub file_format_is_abc{
   my $filename = shift;
   my $ok_line_count = 0;
   open my $fh, "<", "$filename" or die "couldn't open $filename for reading.\n";
   for (1..3) {
      my $line = <$fh>;
      last unless($line =~ /^\S+\s+\S+\s+\S+/);
      $ok_line_count++;
   }
   close $fh;
   return ($ok_line_count == 3);
}

sub file_format_is_iie{
   my $filename = shift;
   my $ok_line_count = 0;
   open my $fh, "<", "$filename" or die "couldn't open $filename for reading.\n";
   my $line = <$fh>;
   $ok_line_count++ if($line =~ /^\S+/);
   $line = <$fh>;
   last unless($line =~ /^\s+\S+\s+\S+/);
   $ok_line_count++;
   close $fh;
   return ($ok_line_count == 2);
}
