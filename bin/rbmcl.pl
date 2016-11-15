#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);
use Getopt::Long;
use Graph;

# clustering based on (approximate) reciprocal best matches
# input is abc file (i.e. each line has id1 id2 e-value)
# make a graph with edges joining reciprocal best matches
# clusters are the biconnected components of the graph.

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

my $weed = 1; # weed out size 2 clusters if both seqs are also in bigger clusters.
my $clusters_filename = 'clusters';
my $abc_filename = undef;
my $gg_filename = '/home/tomfy/Aug2015multispeciesquery/55set.gg';
my $queryspecies_list_filename;
# store best match of each species

GetOptions(
	   'abc_filename=s'           => \$abc_filename, #
	   'gg_filename=s'          => \$gg_filename, # 
           'output_filename=s' => \$clusters_filename,
           'qspecies_filename=s' => \$queryspecies_list_filename,
           'weed!' => \$weed,
);

#my $abc_filename = shift; 
#my $gg_filename = shift || '/home/tomfy/Aug2015multispeciesquery/55set.gg';
my $geneid_sp = store_gg_info("$gg_filename");
print STDERR "\n", "# Done storing gg info.\n";
my %qspecies_count =  ();
#  'Medicago_truncatula' => 0, 'Lotus_japonicus' => 0,
#  'Carica_papaya' => 0, 'Solanum_lycopersicum' => 0,
#  'Oryza_sativa' => 0, 'Phoenix_dactylifera' => 0,
# );
open my $fh_qsp, "<", "$queryspecies_list_filename" or die "Couldn't open $queryspecies_list_filename \n";
while (<$fh_qsp>) {
   next if(/^\s*#/); 
   if (/^\s*(\S+)/) {
      $qspecies_count{$1} = 0;
   }
}

my %id1_id2ev = ();
my %id1__sp_id2 = ();
my %id_selfev = ();
my %id_bccos = (); # key: id, val: array ref of refs to bcc (biconnected component) objs it belongs to. 
my %id_ref = (); # key is an id (string), value is ref to scalar containing that string
my $G = Graph::Undirected->new();

### read in blast match info from abc file:

open my $fh_abc, "<", "$abc_filename" or die "Couldn't open $abc_filename for reading. Exiting. \n";
my $old_id1 = 'not_an_actual_id';
while (my $line = <$fh_abc>) {
   my ($id1, $id2, $ev) = split(" ", $line);
   $id_ref{$id1} = \$id1 if(!exists $id_ref{$id1});
   $id_ref{$id2} = \$id2 if(!exists $id_ref{$id2});
   
   if ($id1 ne $old_id1) {
      $id1_id2ev{$id1} = {};
      $id_bccos{$id1} = [];
      $G->add_vertex($id1);
      for my $sp (keys %qspecies_count) {
         $qspecies_count{$sp} = 0;
      }

   }
   $id1_id2ev{$id1}->{$id2} = $ev;
   my $sp1 = $geneid_sp->{$id1};
   next if(!exists $qspecies_count{$sp1}); # skip if not one of the specified species.
   my $sp2 = $geneid_sp->{$id2};
   next if(!exists $qspecies_count{$sp2}); # skip if not one of the specified species.
 
   if ($qspecies_count{$sp2} == 0) { # skip all but first (best) match.
      my $idref = (exists $id_ref{$id2})? $id_ref{$id2} : die "No ref to $id2 stored in id_ref hash.\n";
      $id1__sp_id2{$id1}->{$sp2} = $idref; # {evalue => $ev, id => $idref}; # "$id2 $ev";
      $qspecies_count{$sp2} = 1;
  #    print STDERR "match: $id1  $id2 \n";
   }
   $old_id1 = $id1;
}
print STDERR "# Done reading in abc data. ";
print STDERR scalar keys %id1__sp_id2, "  query ids. \n";
####### Done reading in abc data #####

################## add edges to the graph, joining reciprocal best matches.

my @qids = sort keys %id1__sp_id2;
for my $qid (@qids) {
   print STDERR "ZZZZQ: ", $qid, "\n";
   my $sp_id2ref = $id1__sp_id2{$qid};
   my $qsp = $geneid_sp->{$qid};
   while (my ($sp, $id2ref) = each %$sp_id2ref) { # get the best matches of each species to $qid
   #   print STDERR "AAAAA: $qsp  $sp \n";
      next if($sp eq $qsp); # don't add edges joining 2 seqs of same species.
   #   print STDERR "ASDFGH\n";
      my $id2 = ${$id2ref};
      print STDERR "DDDDDDDDDDD: $qsp  $qid     $sp  $id2   ", ${$id1__sp_id2{$id2}->{$qsp}}, "\n" if(exists$id1__sp_id2{$id2}->{$qsp});
      next if($id2 eq $qid); # don't add edge connecting node to itself.
      if (exists $id1__sp_id2{$id2}->{$qsp} and ( ${$id1__sp_id2{$id2}->{$qsp}} eq $qid )) { # check for a reciprocal best match
        print  STDERR "adding edge  between $qid  and  $id2 \n";
         $G->add_edge($qid, $id2);
      }
   }
}
undef %id1__sp_id2;
print STDERR "Done adding edges to graph.\n";
print STDERR "Graph has ", scalar $G->vertices(), " vertices, and ", scalar $G->edges(), " edges.\n";
###### Done adding edges to graph ###############

####################### Get biconnected components ########
my @biconnected_components = $G->biconnected_components();
print STDERR "Done getting biconnected components.\n";
######### Done getting biconnected components #############

######### Store array of biconnected-component objects, and %id_bccos hash with the bcc object or objects each id belongs to
my $n_biconn_comps = scalar @biconnected_components;
my @bccobjs = ();
my $sum_cc_sizes = 0;
my $cluster_string = '';
while (my ($index, $ccomponent) = each @biconnected_components) { # loop over connected components (i.e. families)
   # $ccomponent is an array ref of vertices
   $sum_cc_sizes += scalar @$ccomponent;
   my $min_degree = 1000000;
   my $max_degree = -1;
   my @sverts = sort @$ccomponent;
   for my $v (@sverts) {        #(@$ccomponent) {
      $min_degree = min($min_degree, $G->degree($v));
      $max_degree = max($max_degree, $G->degree($v));
   }
   my $bcc_object = { ids => \@sverts,
                      size => scalar @sverts, min_degree => $min_degree, max_degree => $max_degree };
   push @bccobjs, $bcc_object;
   for my $v (@sverts) {
      if (exists $id_bccos{$v}) {
         push @{$id_bccos{$v}}, $bcc_object; # store the indices of the biconn components to which this vert belongs.
      } else {
         die '$id_bccos{$id} doesnt exist for $id = '. "$v\n";
      }
   }
}
# at this point we have the information we want (biconnected components) from the graph, and don't need
# the graph object itself anymore
undef $G; # will this allow the graph object to be successfully garbage-collected?
undef @biconnected_components;

############### Weed out 'bridge' clusters 
my ($bcc_objs_to_keep, $bridges) = ($weed)? weed_bridges(\@bccobjs, \%id_bccos) : (\@bccobjs, []);
print STDERR "After weed_bridges.  ", scalar @$bridges, " found.\n";
############### Done weeding 'bridge' clusters

################# Output ################
open my $clusters_fh, ">", "$clusters_filename" or die "couldn't open $clusters_filename for writing.\n";
for my $bccobj (@$bcc_objs_to_keep) {
   my $idstr = join(",", @{$bccobj->{ids}});
   print $clusters_fh  "$idstr ", $bccobj->{size}, "  ", $bccobj->{min_degree}, "  ", $bccobj->{max_degree}, "\n";
}
my ($sum_bcc_count, $singleton_count) = print_singles(\%id_bccos, $clusters_fh);
close $clusters_fh;

# open my $singles_fh, ">", "$singles_filename" or die "couldn't open $singles_filename for writing.\n";
# my ($sum_bcc_count, $singleton_count) = print_singles(\%id_bccos, $singles_fh);
# close $singles_fh;

print STDERR "# vertices: ", scalar keys %id_bccos, "\n";
print STDERR "# n size>=2 ccs: $n_biconn_comps   sum cc sizes: $sum_cc_sizes \n";
print STDERR "# n bridges weeded: ", scalar @$bridges, "\n";
print STDERR "#  singletons:  $singleton_count \n";
print STDERR "# avg number of clusters a seq belongs to: ", $sum_bcc_count/(scalar keys %id_bccos), "\n";

###############################################################################################################



###############################################################################################################
# subroutines
###############################################################################################################

# remove biconnected components of size 2 if both vertices are also
# members of biconnected components with size > 2
sub weed_bridges{
   my $bcc_objs = shift;
   my $id_bccobjs = shift;
   my @keeps = ();
   my @bridges = ();
   while (my ($index, $bcco) = each @$bcc_objs) { # loop over connected components
      if ($bcco->{size} > 2) {
         push @keeps, $bcco;
      } else {                  # size 2 ccomponent
         my $both_also_in_big = 1;
         for my $v (@{$bcco->{ids}}) {
            my @bccobjs = @{$id_bccobjs->{$v}};
            my $one_big = 0;
            for my $bcco (@bccobjs) { # loop over all biconn components that $v belongs to.
               if ($bcco->{size} > 2) {
                  $one_big = 1;
                  last;
               }
            }
            if (!$one_big) {
               $both_also_in_big = 0;
            }
         }
         if($both_also_in_big){
            push @bridges, $bcco;
         }else{
         push @keeps, $bcco;
      }
      }
   }
   return (\@keeps, \@bridges);
}

sub print_singles{
   my $id_bccos = shift;
   my $fh = shift;
   my $sum_bcc_count = 0;
   while (my ($i, $bccos) = each %$id_bccos) {
      my $bcc_count = scalar @$bccos;
      $sum_bcc_count += ($bcc_count > 0)? $bcc_count : 1;
      if ($bcc_count == 0) {
         print $fh "$i \n";
         $singleton_count++;
      }
   }
   return ($sum_bcc_count, $singleton_count);
}
