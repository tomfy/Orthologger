package PipeSegRbmcl;
use strict;
use parent 'PipelineSegment';
use List::Util qw (min max sum);
use Cwd qw (getcwd abs_path);
use Pod::Usage;
use Graph;

# clustering based on (approximate) reciprocal best matches
# input is abc file (i.e. each line has id1 id2 e-value)
# make a graph with edges joining reciprocal best matches
# clusters are the biconnected components of the graph.

sub new{
   my $class = shift;
   my $segment_name = shift;
   my $pl_obj = shift;   # the pipeline object this segment belongs to
   print STDERR "Top of PipeSegRbmcl constructor.\n";
   my $self = $class->SUPER::new(); # PipelineSegment obj.
   $self->init(segment_name => $segment_name,
               pipeline => $pl_obj,
               weed => 1,
              );
   $self->hide('id1_id2ev');
   return $self;
}

sub run{
   my $self = shift;

   my $pl_obj = $self->get('pipeline');
   my $fh_progress = $pl_obj->get('fh_progress');

   print STDERR "top of PipeSegRbmcl->run; ", getcwd(), " ", abs_path(getcwd()), "\n";

   my $out_dir;

   print STDERR "In PSRbmcl::run predecessor: ", ($self->get('predecessor')) // 'undef', "\n"; #exit;
   if (defined $self->get('predecessor')) {
      print STDERR "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA predecessor:  ", $self->get('predecessor'), "\n";
      $out_dir = $self->get('output_dir');
      print STDERR "output dir: $out_dir \n";
      mkdir $out_dir unless(-d $out_dir); # create dir if doesn't exist
      chdir($out_dir);
      print STDERR "pwd: ", abs_path('./'), "\n";
   } else { # should be in dir (created by hand) which is subdir of blast output dir. should have *.m8 files
      $out_dir = getcwd();
   }
exit;
   # ***   make iie files   ***
   my $m8_filename_str = `ls ../*.m8`;
   die "ls ../*.m8 finds no files.\n" if($m8_filename_str =~ /No such file or directory/);
   print STDERR "ZZZ  m8_filename_str: $m8_filename_str \n";
   my @m8_filenames = split(" ", $m8_filename_str);
   my $id1__id2_ev = $self->m8toiie(\@m8_filenames);

   # ***   write .iie file (needed?)   ***
   my $id1id2ev_string = '';
   for my $id1 ($id1__id2_ev->keys()) {
      my $id2ev = $id1__id2_ev->get($id1);
      $id1id2ev_string .= "$id1\n";
      for my $id2 ($id2ev->keys) {
         $id1id2ev_string .= "  $id2  " . $id2ev->get($id2) . "\n";
      }
   }
   open my $fhiie, ">", 'matches.iie' or die "Couldn't open iie for writing.\n";
   print $fhiie $id1id2ev_string;
   close $fhiie;
   $self->set('matches_filename' => 'matches.iie');

   # ***   rbm clustering   ***
   $self->rbmcl();

   # ***   wrap up   ***
   $self->set(state => 'complete');
 #  $self->set_successor_state('next');
}


sub m8toiie{ # from arrayref of m8 format filenames to Hash::Ordered of id1/Hash::Ordered pairs.
   my $self = shift;
   my $filenames = shift;       # array ref of .m8 filenames
   my $id_sp = shift || undef; # optional hashref giving id/sp association; if defined, skip id1's which are not keys.
   my $qsp = shift || undef; # optional hashref whose keys are the query species, skips any id1 not a key of $qsp.

   my $id1__id2_ev = Hash::Ordered->new(); # ();        #
   for my $fname (@$filenames) {
      my ($id1, $id2, $evalue);
      open my $fh, "<", "$fname" or die "Couldn't open file $fname for reading.\n";
      while (my $line = <$fh>) {
         next if($line =~ /^\s*#/); # skip comments
         my @cols = split(' ', $line);
         next if(scalar @cols != 12); # there should be 12 columns
         next if( (defined $id1 and defined $id2)  and  ($cols[0] eq $id1)  and  ($cols[1] eq $id2) );
         ($id1, $id2, $evalue) = @cols[0,1,10];
         next if((defined $id_sp and defined $qsp)  and (!exists $id_sp->{$id1} or !exists $qsp->{$id_sp->{$id1}}));
         if (!defined $id1__id2_ev->get($id1)) {
            $id1__id2_ev->set($id1 => Hash::Ordered->new());
         }
         $id1__id2_ev->get($id1)->set($id2 => $evalue);
      }
      close $fh;
   }
   $self->set(id1_id2ev => $id1__id2_ev);
   return $id1__id2_ev;
}

sub rbmcl{
   my $self = shift;
   my $iie_filename = $self->get('matches_filename');

   my $pl_obj = $self->get('pipeline');
   my %qspecies_count =  ();
   for ($pl_obj->get('querytaxon_inputpath')->keys()) {
      $qspecies_count{$_} = 0;
   }

   while (my($k, $v) = each %qspecies_count) {
      print STDERR "ABCDEF: $k $v \n";
   }
   if (!defined $self->get('clusters_filename')) {
      my $clusters_filename = $iie_filename;
      $clusters_filename =~ s/([.]iie)/.qclusters/; # X.iie -> X.qclusters
      $self->set('clusters_filename' => $clusters_filename);
   }
   my ($id1__sp_id2, $id_bccos) = 
     #   $self->read_iie($iie_filename);
     $self->iie2x();
   ################## add edges to the graph, joining reciprocal best matches.
   my $G = Graph::Undirected->new();
   my @qids = sort keys %$id1__sp_id2;
   for my $qid (@qids) {
      $G->add_vertex($qid);
   }
   for my $qid (@qids) {
      my $sp_id2 = $id1__sp_id2->{$qid};
      my $qsp = $pl_obj->get('seqid_species')->{$qid};
#      print STDERR "QSP QID: $qsp  $qid.\n";
      while (my ($sp, $id2) = each %$sp_id2) { # get the best matches of each species to $qid
#         print STDERR "QSP TSP TID: $qsp  $sp $id2 \n";
         next if($sp eq $qsp); # don't add edges joining 2 seqs of same species.
         if (exists $id1__sp_id2->{$id2}->{$qsp}) {
 #           print STDERR '$id1__sp_id2->{$id2}->{$qsp}: ', $id1__sp_id2->{$id2}->{$qsp}, " $qid \n";
            if ($id1__sp_id2->{$id2}->{$qsp} eq $qid ) { # check for a reciprocal best match
               $G->add_edge($qid, $id2);
 #              print STDERR "add edge between: $qid, $id2!!!!\n";
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
      for my $v (@sverts) {     #(@$ccomponent) {
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
   my ($bcc_objs_to_keep, $bridges) = ($self->get('weed'))? weed_bridges(\@bccobjs, $id_bccos) : (\@bccobjs, []);
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

   my $clusters_filename = $self->get('clusters_filename');
   open my $clusters_fh, ">", "$clusters_filename" or die "couldn't open $clusters_filename for writing.\n";
   print $clusters_fh join("", sort @out_lines);
   close $clusters_fh;

   print STDERR "# vertices: ", scalar keys %$id_bccos, "\n";
   print STDERR "# n size>=2 ccs: $n_biconn_comps   sum cc sizes: $sum_cc_sizes \n";
   print STDERR "# n bridges weeded: ", scalar @$bridges, "\n";
   print STDERR "#  singletons:  $singleton_count \n";

   ###############################################################################################################
}

###############################################################################################################
# subroutines (non-methods)
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
   my $singleton_count = 0;
   while (my ($i, $bccos) = each %$id_bccos) {
      my $bcc_count = scalar @$bccos;
      if ($bcc_count == 0) {
         push @singles, sprintf("%s\n", $i);
         $singleton_count++;
      }
   }
   return (\@singles, $singleton_count);
}

sub read_iie{
   my $self = shift;
   my $iie_filename = shift;

   my $id1__sp_id2 = {};
   my $id_bccos = {}; # key: id, val: array ref of refs to bcc (biconnected component) objs it belongs to. 
   my $geneid_sp = $self->get('pipeline')->get('seqid_species');
   if (!defined $geneid_sp) {
      $geneid_sp = $self->get('pipeline')->get('segments')->get('gg')->get('seqid_species');
      print STDERR "Number of gene ids in geneid_sp: ", (defined $geneid_sp)? scalar keys %$geneid_sp : 'undef', "\n";
   }
   ### read in blast match info from iie file:
   $iie_filename = abs_path($iie_filename);
   print STDERR "$iie_filename \n";
   open my $fh_iie, "<", "$iie_filename" or die "Couldn't open $iie_filename for reading. Exiting. \n";
   my ($id1, $id2, $ev);
   #   my @queryspecies = $self->get('pipeline')->get('querytaxon_inputpath')->keys();
   my $qtip_obj = $self->get('pipeline')->get('querytaxon_inputpath');
   my %qspecies_count = ();     # map { $_ => 0 } @queryspecies;
   while (my $line = <$fh_iie>) {
      if ($line =~ /^(\S+)/) {  # this is a line with just a query id.
         $id1 = $1;
         $id_bccos->{$id1} = [];
         %qspecies_count = map { $_ => 0 } $qtip_obj->keys();
 #        print STDERR "Query taxa: ", join("; ", keys %qspecies_count), "\n"; #exit;
      } elsif ($line =~ /^\s+(\S+)/) { # \s+(\S+)/) {
         $id2 = $1;
         my $sp1 = $geneid_sp->{$id1};
         next if(!exists $qspecies_count{$sp1}); # skip if not one of the specified species.
         my $sp2 = $geneid_sp->{$id2};
  #       print STDERR "sp1 sp2: $sp1 $sp2 \n";
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

sub iie2x{
   my $self = shift;
   my $id1_id2ev = $self->get('id1_id2ev'); # Hash::Ordered of id1/Hash::Ordered pairs.

   my $id1__sp_id2 = {};
   my $id_bccos = {}; # key: id, val: array ref of refs to bcc (biconnected component) objs it belongs to.
   my $geneid_sp = $self->get('pipeline')->get('seqid_species');
   my @querytaxa = $self->get('pipeline')->get('querytaxon_inputpath')->keys();
   my %qspecies_count = ();
   #for(
   my ($id1, $id2, $ev);
   my $iter = $id1_id2ev->iterator();
   while (my ($id1, $id2ev) = $iter->()) { # each %{$id1_id2ev}){
      $id_bccos->{$id1} = [];
      %qspecies_count = map { $_ => 0 } @querytaxa; # $qtip_obj->keys();

  #    print STDERR "Query taxa: ", join("; ", keys %qspecies_count), "\n"; #exit;
 #     print STDERR "id1: $id1 \n";
      for my $id2 ($id2ev->keys()) {
   #      print STDERR "   $id2 \n";
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

1;
