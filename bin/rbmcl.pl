#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);
use Getopt::Long;
use Pod::Usage;
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
use TomfyMisc qw 'file_format_is_abc_iie run_quicktree run_fasttree run_phyml store_gg_info timestring ';

my $weed = 1; # weed out size 2 clusters if both seqs are also in bigger clusters.
my $clusters_filename = undef;
my $matches_filename = undef;
my $gg_filename = '/home/tomfy/Aug2015multispeciesquery/55set.gg';
my $queryspecies_list_filename;
my $help = undef;
# store best match of each species

GetOptions(
           'help|?:i' => \$help, # sub {my ($optname, $optval) = @_, return ($optval == 0}, # \$help,
	   'matches_filename=s'           => \$matches_filename, #
           #     'iie_filename=s' => \$iie_filename,
	   'gg_filename=s'          => \$gg_filename, #
           'output_filename=s' => \$clusters_filename,
           'qspecies_filename=s' => \$queryspecies_list_filename,
           'weed!' => \$weed,
          );

my $geneid_sp = store_gg_info("$gg_filename");
print STDERR "\n", "# Done storing gg info.\n";

my %qspecies_count =  ();
open my $fh_qsp, "<", "$queryspecies_list_filename" or die "Couldn't open $queryspecies_list_filename \n";
while (<$fh_qsp>) {
   next if(/^\s*#/); 
   if (/^\s*(\S+)/) {
      $qspecies_count{$1} = 0;
   }
}

while(my($k, $v) = each %qspecies_count){
   print STDERR "ABCDEF: $k $v \n";
}
if (!defined $clusters_filename) {
   $clusters_filename = $matches_filename;
   $clusters_filename =~ s/([.]iie|[.]abc)/.qclusters/; # X.iie -> X.qclusters
}
my ($id1__sp_id2, $id_bccos);
my $matches_format = file_format_is_abc_iie($matches_filename);
if ($matches_format eq 'iie') {
   ($id1__sp_id2, $id_bccos) = read_iie($matches_filename);
} elsif ($matches_format eq 'abc') {
   ($id1__sp_id2, $id_bccos) = read_abc($matches_filename);
} else {
   die "matches file $matches_filename has unrecognized format.\n";
}
################## add edges to the graph, joining reciprocal best matches.
my $G = Graph::Undirected->new();
my @qids = sort keys %$id1__sp_id2;
for my $qid (@qids) {
   $G->add_vertex($qid);
}
for my $qid (@qids) {
   my $sp_id2 = $id1__sp_id2->{$qid};
   my $qsp = $geneid_sp->{$qid};
   print STDERR "$qsp  $qid.\n";
   while (my ($sp, $id2) = each %$sp_id2) { # get the best matches of each species to $qid
      print STDERR "$qsp  $sp $id2 \n";
      next if($sp eq $qsp); # don't add edges joining 2 seqs of same species.
      if (exists $id1__sp_id2->{$id2}->{$qsp}) {
         if ($id1__sp_id2->{$id2}->{$qsp} eq $qid ) { # check for a reciprocal best match
            $G->add_edge($qid, $id2);
            print STDERR "add edge between: $qid, $id2!!!!\n";
         }
      }
   }
}
undef %$id1__sp_id2;
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
      if (exists $id_bccos->{$v}) {
         push @{$id_bccos->{$v}}, $bcc_object; # store the indices of the biconn components to which this vert belongs.
      } else {
         die '$id_bccos->{$id} doesnt exist for $id = '. "$v\n";
      }
   }
}
# at this point we have the information we want (biconnected components) from the graph, and don't need
# the graph object itself anymore
undef $G; # will this allow the graph object to be successfully garbage-collected?
undef @biconnected_components;
############### Weed out 'bridge' clusters 
my ($bcc_objs_to_keep, $bridges) = ($weed)? weed_bridges(\@bccobjs, $id_bccos) : (\@bccobjs, []);
print STDERR "After weed_bridges.  ", scalar @$bridges, " found.\n";
############### Done weeding 'bridge' clusters

################# Output ################
my @out_lines = ();
for my $bccobj (@$bcc_objs_to_keep) {
   my $idstr = join(",", @{$bccobj->{ids}});
   push @out_lines, "$idstr ". $bccobj->{size} . "  " . $bccobj->{min_degree} . "  " . $bccobj->{max_degree} . "\n";
}
my ($singles_aref, $singleton_count) = singles_string($id_bccos);
push @out_lines, @{$singles_aref};

open my $clusters_fh, ">", "$clusters_filename" or die "couldn't open $clusters_filename for writing.\n";
print $clusters_fh join("", sort @out_lines);
close $clusters_fh;

print STDERR "# vertices: ", scalar keys %$id_bccos, "\n";
print STDERR "# n size>=2 ccs: $n_biconn_comps   sum cc sizes: $sum_cc_sizes \n";
print STDERR "# n bridges weeded: ", scalar @$bridges, "\n";
print STDERR "#  singletons:  $singleton_count \n";

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
         if ($both_also_in_big) {
            push @bridges, $bcco;
         } else {
            push @keeps, $bcco;
         }
      }
   }
   return (\@keeps, \@bridges);
}

sub singles_string{
   my $id_bccos = shift;
   my @singles = ();
   while (my ($i, $bccos) = each %$id_bccos) {
      my $bcc_count = scalar @$bccos;
      if ($bcc_count == 0) {
         push @singles, sprintf("%s\n", $i);
         $singleton_count++;
      }
   }
   return (\@singles, $singleton_count);
}

sub read_abc{
   my $abc_filename = shift;

   my $id1__sp_id2 = {};
   my $id_bccos = {}; # key: id, val: array ref of refs to bcc (biconnected component) objs it belongs to. 

   ### read in blast match info from abc file:

   open my $fh_abc, "<", "$abc_filename" or die "Couldn't open $abc_filename for reading. Exiting. \n";
   my $old_id1 = 'not_an_actual_id';
   while (my $line = <$fh_abc>) {
      my ($id1, $id2, $ev) = split(" ", $line);
      if ($id1 ne $old_id1) {
         $id_bccos->{$id1} = [];
         for my $sp (keys %qspecies_count) {
            $qspecies_count{$sp} = 0;
         }
      }

      my $sp1 = $geneid_sp->{$id1};
      next if(!exists $qspecies_count{$sp1}); # skip if not one of the specified species.
      my $sp2 = $geneid_sp->{$id2};
      next if(!exists $qspecies_count{$sp2}); # skip if not one of the specified species.

      if ($qspecies_count{$sp2} == 0) { # skip all but first (best) match.
         $id1__sp_id2->{$id1}->{$sp2} = $id2; # $idref; # {evalue => $ev, id => $idref}; # "$id2 $ev";
         $qspecies_count{$sp2} = 1;
      }
      $old_id1 = $id1;
   }
   print STDERR "# Done reading in abc data. ";
   print STDERR scalar keys %$id1__sp_id2, "  query ids. \n";
   ####### Done reading in abc data #####
   return ($id1__sp_id2, $id_bccos);
}


sub read_iie{
   my $iie_filename = shift;

   my $id1__sp_id2 = {};
   my $id_bccos = {}; # key: id, val: array ref of refs to bcc (biconnected component) objs it belongs to. 

   ### read in blast match info from iie file:
   print STDERR "$iie_filename \n";
   open my $fh_iie, "<", "$iie_filename" or die "Couldn't open $iie_filename for reading. Exiting. \n";
   my ($id1, $id2, $ev);
   while (my $line = <$fh_iie>) {
      if ($line =~ /^(\S+)/) {
         $id1 = $1;
         $id_bccos->{$id1} = [];
         for my $sp (keys %qspecies_count) {
            $qspecies_count{$sp} = 0;
         }
      } elsif ($line =~ /^\s+(\S+)/) { # \s+(\S+)/) {
         $id2 = $1;
         my $sp1 = $geneid_sp->{$id1};
         next if(!exists $qspecies_count{$sp1}); # skip if not one of the specified species.
         my $sp2 = $geneid_sp->{$id2};
         next if(!exists $qspecies_count{$sp2}); # skip if not one of the specified species.

         if ($qspecies_count{$sp2} == 0) { # skip all but first (best) match.
            $id1__sp_id2->{$id1}->{$sp2} = $id2; # {evalue => $ev, id => $idref}; # "$id2 $ev";
            $qspecies_count{$sp2} = 1;
         }
      }
   }
   print STDERR "# Done reading in iie data. ";
   print STDERR scalar keys %$id1__sp_id2, "  query ids. \n";
   ####### Done reading in iie data #####
   return ($id1__sp_id2, $id_bccos);
}


__END__


=head1 NAME

    rbmcl.pl - cluster query ids, by constructing graph with reciprocal best matches joined
                  by edges. clusters are biconnected components of graph.

=head1 SYNOPSIS

    rbmcl.pl  -gg_filename <filename> -matches_filename <filename> -qspecies_filename <filename> [options]
     Options:
       -gg_filename  gene-genome association file; specifies genes in each genome.
       -matches_filename  blast results in abc or iie format. Required - no default.
       -qspecies_filename  file specifying the query species, 1 species per line. Required - no default.
       -output_filename Default: construct from input filename, with form *_fams.iis
       -weed  maximum number of matches to include. Default: 400

=head1 OPTIONS

=over 2

=item B<-gg_filename> 

             gene-genome association file. 1 genome per line; e.g. 
              'Amborella_trichopoda: evm_27.model.AmTr_v1.0_scaffold00001.498 evm_27.model.AmTr_v1.0_scaffold00001.491 (etc.)' 

=item B<-matches_filename>

    File with blast output information. Format is either abc (id1 id2 e-value on each line), or
     iie (id1 on first line of family, then ' id2 e-value' (with initial space) on each line
     for matches to id1.)

=item B<-qspecies_filename>

    File specifying which species are query species. 1 species name per line (e.g. Arabidopsis_thaliana)

=item B<-output_filename>

    Name to use for output file. Default is to truncate '.abc' or '.iie' from input filename,
     and append '.qclusters'.

=item B<-max_fam_size>

    -weed or -noweed . Default is weed, meaning that size 2 clusters with both members also belonging
              to larger clusters are eliminated. 

=back

=head1 DESCRIPTION

B<rbmcl.pl> Given matches from blast, and a list of query species, looks for reciprocal best matches among
  the query sequences, constructs a graph whose nodes are query sequences, and with edges joining pairs of 
  nodes which are reciprocal best matches, and then finds the biconnected components of the graph. A biconnected
  component is a set of nodes which is connected, and which remains connected if any ONE of its nodes is deleted.
  (I.e. the remaining nodes are connected.)
  A node can belong to > 1 biconnected component. E.g. imagine two triangles which share a vertex. Each triangle
  is a biconnected component, and the shared vertex belongs to both; the whole (5 vertices) is connected, but 
  not biconnected, because the deletion of the shared vertex results in two separate connected components.


=cut

