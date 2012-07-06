#!/usr/bin/perl 

use strict;
use Getopt::Std;
use List::Util qw ( min max sum );

use lib '/home/tomfy/Orthologger/lib';
use Overlap;
use Orthologger;
use TlyUtil qw ( order_newick );
use Mrbayes;

#use lib '/home/tomfy/cxgn/cxgn-corelibs/lib';
use CXGN::Phylo::File;
use CXGN::Phylo::Parser;

# use Devel::Cycle; # for finding circular refs - cause of memory leaks.
# find_cycle($test); # to find circular refs in $test (which might be an object, e.g.)


# find orthologs, with levels of support - either bootstrap support, or
# bayesian posterior probability from MrBayes.

# read in an alignment file. get an overlap object.
# repeat $n_bootstrap times:
# 1) get bootstrap overlap
# 2) get tree (using fasttree, quicktree (default) or clearcut)/
# 3) optionally reroot this tree in various ways (midpoint, mindl, minvar ...)
# 4) get orthologs for this tree
# keep track of how many times a pair is predicted orthologous
# report bootstrap support for each ortholog pair
# omitting those below threshold

my $seed_increment = 1000; # if seed is given on c.l. increment by this for each bootstrap.
my $do_set_error = 0; # 0 speeds up parsing by skipping many calls to set_error.
my $clearcut_spacer = ' '; # this goes between the '>' and the id in fasta
# input to clearcut. Otherwise clearcut sometimes puts funny characters in its newick output.
my $analyze_mr_bayes = 0; # if true, look for MrBayes output to analyze.

print "# $0 ", join(" ", @ARGV), "\n";

use vars qw($opt_i $opt_s $opt_S $opt_f $opt_T $opt_N $opt_q $opt_Q $opt_k $opt_n $opt_r $opt_m);
# -i <input_alignment_filename>     (name of input alignment file, fasta format)
# -s <species_tree_filename>    (file containing species tree in newick format.)
# -S <seed>  (rng seed, default: get from clock)
# -f <fraction> (fraction of non-gap chars required in alignment column to keep it.)
# -T <tree finding method> (NJ, or ML. Use NJ, or ML to infer tree. Default: NJ)
# -N <n>    (number of bootstrap data sets to create and process. Default is 1.)
# -q <species>  (query species. Default is show all.)
# -Q <sequence_id regex> (query id. Show sequences whose id matches this as a regular expression.)
# -k    (kimura correction to distances. default is no correction.)
# -n <n>  (produce this many trees for the actual data set. default is 1. ignored for NJ, ML)
# -r <reroot type> (mindl midpoint minvar, none, etc. Default is mindl)
# -m <min bootstrap support> (Output only orthologs having >= this bootstrap supportor posterior probability; if > 1 interpret as %)

# typical usage: perl bootstrap_ortholog.pl -i $align_filename -T ML -k -n 1 -N 100 -S 12345 -r mindl -m 0.015 -q castorbean > outfile

# get options
getopts("i:s:S:f:T:N:q:Q:kn:r:m:");

die "usage: bootstrap_ortholog.pl -i <input_filename> [-s species_tree -S <rng seed> -f <overlap fraction> -T <tree infer type>  -N <#bootstrap data sets> -q <query sequence id> -k -n <# rng trees from each data set> i" unless($opt_i and -f $opt_i);

my $input_file = $opt_i;
my $output_tree_filename = $input_file;
$output_tree_filename =~ s/fasta/newick/;
$output_tree_filename .= '.newick' if(! $output_tree_filename =~ /newick/);
my $output_orthologs_filename = $input_file;
$output_orthologs_filename =~ s/fasta/orthologs/;
$output_orthologs_filename .= '.orthologs' if(! $output_orthologs_filename =~ /orthologs/);



my $n_bootstrap = $opt_N || 0; # default is just do the actual data, not bootstraps.
my ($n_NJ_bootstrap, $n_ML_bootstrap) = (0, 0);
if($n_bootstrap =~ /(\d+),(\d+)/){
  $n_NJ_bootstrap = $1;
$n_ML_bootstrap = $2;
}elsif($n_NJ_bootstrap =~ /(\d+)/){
  $n_NJ_bootstrap = $n_bootstrap;
  $n_ML_bootstrap = $n_bootstrap;
}else{
  warn "opt_N: $opt_N; invalid. No bootstraps will be performed.\n";
  $n_bootstrap = 0;
}

my $nongap_fraction = $opt_f || 0.8;

my $query_species = ($opt_q)? $opt_q: undef;
my $query_id_regex = ($opt_Q)? $opt_Q: undef;

my $bootstrap_seed = ($opt_S)? $opt_S: undef;
my $tree_find_seed = ($opt_S)? $opt_S + $seed_increment: undef;

my $clearcut_base_cl = 'clearcut --stdin --stdout -a -P ';
#$clearcut_base_cl = 'MrRogers --stdin --stdout -a -P ';
my $fasttree_base_cl = 'FastTree -wag -gamma -bionj -quiet '; # -noml -nosupport';

my $n_rep = ($opt_n)? $opt_n : 1;
my $tree_find_method = (defined $opt_T)? $opt_T : 'NJ';
if ($tree_find_method eq 'NJ') {
  $clearcut_base_cl .= '-N ';
  $n_rep = 1;
}
my $min_bs_support = (defined $opt_m)? $opt_m : 0.0; #default is keep all.
if ($min_bs_support > 1) {	# interpret as percentage if > 1.
  $min_bs_support *= 0.01;
}
$clearcut_base_cl .= ($opt_k)? '-k ': '';
my $quicktree_distance_correction = ($opt_k)? 'kimura' : 'none';

#### Get the alignment:
my $align_string =  `cat $input_file`;
# fixes to $align_string:
$align_string =~ s/IMGA[|]/IMGA_/g; #pipes in id cause problem; replace '|' with '_'.
#$align_string =~ s/(>[^|]+)[|][^\n]+\n/$1\n/; # delete from first pipe to end of line.
$align_string =~ s/(>[^|]+)[|][^\n]*\n/$1\n/g; # delete from first pipe to end of line.

my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this.
$align_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg; # to make clearcut happy.
# print "align string: $align_string\n";


# construct an overlap object.
my $overlap_obj = Overlap->new($align_string, $nongap_fraction, $bootstrap_seed);

### Setup for orthologger.
# set up cl for orthologger. options $opt_s, $opt_r are relevant.
my $reroot_method = 'mindl';
if (defined $opt_r) {
  if ($opt_r eq 'none') {
    $reroot_method = '';
  } elsif ($opt_r eq 'mindl' or
	   $opt_r eq 'midpoint' or
	   $opt_r eq 'minvar' or
	   $opt_r eq 'maxmin') {
    $reroot_method = "$opt_r";
  } else {
    warn "reroot option $opt_r unknown. using default: mindl\n";
  }
}
# default species tree: 13-species tree:
my $species_newick = "(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
# get the species tree if specified on c.l.
my $species_tree_file_name;
if (defined $opt_s) {
  if (-f $opt_s) {		# read species tree newick from file
    my $species_tree_file = CXGN::Phylo::File->new($opt_s);
    $species_newick = $species_tree_file->get_tree_string();
  } else {
    die "species tree file: [$opt_s]; no such file.\n"; 
  }
}
# my $species_tree = CXGN::Phylo::Parse_newick->new($species_newick, $do_set_error)->parse();
print "YYY in ortholog_support.pl about to call Parse_newick. [$do_set_error]\n";
my $sparser = CXGN::Phylo::Parse_newick->new($species_newick, $do_set_error);
my $species_tree = $sparser->parse();
#find_cycle($sparser);
if (!$species_tree) {
  die"Species tree. Parse_newick->parse() failed to return a tree object. Newick string: "
    . $species_newick."\n";
}
print "# species tree: \n# $species_newick \n";
# find_cycle($species_tree);

my %idpair_orthocount_all = (); # includes results for both actual and bootstrap data.
my %idpair_actual_orthocountNJ = (); # results for actual data only.
my %idpair_actual_orthocountML = (); # results for actual data only.
my %idpair_bs_orthocountNJ = (); # results for ML analyzed bootstrap data.
my %idpair_bs_orthocountML = (); # results for ML analyzed bootstrap data.
my %idpair_orthocountMB = (); # Posterior prob. results from MrBayes.

#my $spacer = $clearcut_spacer; #($tree_find_method eq 'ML')? '' : $clearcut_spacer;
# Analyze actual data:
my $overlap_fasta_string = $overlap_obj->overlap_fasta_string('');
my $clearcut_overlap_fasta_string = $overlap_obj->overlap_fasta_string($clearcut_spacer);
my $overlap_length = $overlap_obj->get_overlap_length();
my $clearcut_cl = $clearcut_base_cl;
my $fasttree_cl = $fasttree_base_cl;
$clearcut_cl .= ($opt_S)? " -s $tree_find_seed ": '';
$fasttree_cl .= ($opt_S)? " -seed $tree_find_seed ": ''; 
#$clearcut_cl .= ($opt_n)? " -n $n_rep ": ''; # actual data, do multiple runs (RNJ or ML)
print "# overlap length: $overlap_length.\n";
#print "# clearcut command line: $clearcut_cl \n";
print "# fasttree command line: $fasttree_cl \n";

my $newicks_out = # run_clearcut($clearcut_overlap_fasta_string, $clearcut_cl);
  run_quicktree($overlap_fasta_string, $quicktree_distance_correction);

my $first = 1;
my $rerooted_gene_tree_newick;
foreach my $gene_tree_newick ( split("\n", $newicks_out) ) {
  next if($gene_tree_newick =~ /^NJ/);
  $gene_tree_newick =~ s/;\s*$//; # this is gene tree newick on one line.
 
  my $orthologger_obj = Orthologger->new({'gene_tree_newick' => $gene_tree_newick, 'species_tree' => $species_tree, 'reroot_method' => $reroot_method, 'query_species' => $query_species, 'query_id_regex' => $query_id_regex});
  print "# Orthologger object constructed.\n";
  if ($first) {
    $rerooted_gene_tree_newick = $orthologger_obj->get_gene_tree()->generate_newick();
    print "# Actual data rerooted NJ gene tree newick:\n# $rerooted_gene_tree_newick\n";
    $first = 0;
  }
  my $orthologger_outstring = $orthologger_obj->ortholog_result_string();

  store_orthologger_out($orthologger_outstring, [\%idpair_orthocount_all, \%idpair_actual_orthocountNJ]);
  $orthologger_obj->decircularize();
  #find_cycle($orthologger_obj);
}

if ($tree_find_method eq 'ML') {
  my $gene_tree_newick = run_fasttree($overlap_fasta_string, $fasttree_cl);
  my $first = 1;
  my $rerooted_gene_tree_newick;
  #foreach my $gene_tree_newick ( split("\n", $newicks_out) ) {
  #  next if($gene_tree_newick =~ /^NJ/);
  $gene_tree_newick =~ s/;\s*$//; # this is gene tree newick on one line.
  my $orthologger_obj = Orthologger->new({'gene_tree_newick' => $gene_tree_newick, 'species_tree' => $species_tree, 'reroot_method' => $reroot_method, 'query_species' => $query_species, 'query_id_regex' => $query_id_regex});
  if ($first) {
    $rerooted_gene_tree_newick = $orthologger_obj->get_gene_tree()->generate_newick();
    print "# Actual data rerooted $tree_find_method gene tree newick:\n# $rerooted_gene_tree_newick\n";
    print "output tree filename: $output_tree_filename \n";

    open my $fh_tree, ">$output_tree_filename";
    print $fh_tree "# Actual data rerooted $tree_find_method gene tree newick:\n $rerooted_gene_tree_newick\n";
close $fh_tree;

    $first = 0;
  }
  my $orthologger_outstring = $orthologger_obj->ortholog_result_string();

  store_orthologger_out($orthologger_outstring, [\%idpair_orthocount_all, \%idpair_actual_orthocountML]);
  $orthologger_obj->decircularize();
}
print STDERR "Actual data done.\n";
# print STDERR "XXX: ", $rerooted_gene_tree_newick, "\n";

$tree_find_seed += $seed_increment;
# Done with actual data.

# Now the bootstrap data.
for my $i_bs (1..$n_NJ_bootstrap) {
  my $overlap_fasta_string = $overlap_obj->bootstrap_overlap_fasta_string('');
  my $clearcut_overlap_fasta_string = $overlap_fasta_string;
  $clearcut_overlap_fasta_string =~ s/>(\S)/> $1/g; # put in the space
  my $clearcut_cl = $clearcut_base_cl . (($opt_S)? " -s $tree_find_seed " : '');
  my $fasttree_cl = $fasttree_base_cl . ''; # " -noml -nosupport "; 
    $fasttree_cl = $fasttree_cl . (($opt_S)? " -seed $tree_find_seed " : '');

  #print STDERR "$i_bs  $overlap_fasta_string \n";
  #print STDERR "$i_bs  $clearcut_cl \n";

  #  my $clearcut_newicks_out = run_clearcut($overlap_fasta_string, $clearcut_cl);

  my $newicks_out = #run_clearcut($clearcut_overlap_fasta_string, $clearcut_cl);
    run_quicktree($overlap_fasta_string, $quicktree_distance_correction);
  my @ccnewicks = split("\n", $newicks_out);
  foreach my $gene_tree_newick (@ccnewicks) {
    next if($gene_tree_newick =~ /^NJ/);
    $gene_tree_newick =~ s/;\s*$//; # this is gene tree newick on one line.

    my $orthologger_obj = Orthologger->new({'gene_tree_newick' => $gene_tree_newick, 'species_tree' => $species_tree, 'reroot_method' => $reroot_method, 'query_species' => $query_species, 'query_id_regex' => $query_id_regex});

    my $orthologger_outstring = $orthologger_obj->ortholog_result_string();
    $orthologger_obj->decircularize();
    #find_cycle($orthologger_obj);
    store_orthologger_out($orthologger_outstring, [\%idpair_orthocount_all, \%idpair_bs_orthocountNJ]);
  }
  print STDERR "\r                          \r";
  print STDERR "NJ Bootstraps 1..$i_bs done.";
} # end of NJ bootstraps
print "\n";


  if ($tree_find_method eq 'ML') { # also construct/analyze ML trees for the bootstrap
for my $i_bs (1..$n_ML_bootstrap) {
 my $overlap_fasta_string = $overlap_obj->bootstrap_overlap_fasta_string('');
 my $fasttree_cl = $fasttree_base_cl . (($opt_S)? " -seed $tree_find_seed " : '');
    my $gene_tree_newick = run_fasttree($overlap_fasta_string, $fasttree_cl);
    $gene_tree_newick =~ s/;\s*$//; # this is gene tree newick on one line.

    my $orthologger_obj = Orthologger->new({'gene_tree_newick' => $gene_tree_newick, 'species_tree' => $species_tree, 'reroot_method' => $reroot_method, 'query_species' => $query_species, 'query_id_regex' => $query_id_regex});

    my $orthologger_outstring = $orthologger_obj->ortholog_result_string();
    $orthologger_obj->decircularize();
    #find_cycle($orthologger_obj);
    store_orthologger_out($orthologger_outstring, [\%idpair_orthocount_all, \%idpair_bs_orthocountML]);


  $tree_find_seed += $seed_increment;
  print STDERR "\r                          \r";
  print STDERR "ML Bootstraps 1..$i_bs done.";
}
}
print STDERR "\n";		# loop over bootstraps

# do Bayesian analysis
my $alignment_nex_filename = $input_file . ".nex";
$alignment_nex_filename =~ s/[.]fasta//;
#print "$alignment_nex_filename \n";

my $MB_tree_count = process_MB_out($alignment_nex_filename);

my $sum_MB_count = 0;
if($analyze_mr_bayes){ # do Bayesian analysis (using precomputed MrBayes output, if present), if not precomputed skip.
foreach my $gene_tree_newick (keys %$MB_tree_count){
my $MB_count = $MB_tree_count->{$gene_tree_newick};
$sum_MB_count += $MB_count;
# print "XXXXXXX: $MB_count  ", $gene_tree_newick, "\n";
my $d = 1; # distance to use
$gene_tree_newick =~ s/(\S)([,)])/$1:$d$2/g; # put in distances.
$gene_tree_newick =~ s/([)])([,)])/$1:$d$2/g;

 #print "XXXXX [$MB_count]  $gene_tree_newick \n";
 my $orthologger_obj = Orthologger->new({'gene_tree_newick' => $gene_tree_newick, 'species_tree' => $species_tree, 'reroot_method' => $reroot_method, 'query_species' => $query_species, 'query_id_regex' => $query_id_regex});
  # if ($first) {
  #   $rerooted_gene_tree_newick = $orthologger_obj->get_gene_tree()->generate_newick();
  #   print "# Actual data rerooted NJ gene tree newick:\n# $rerooted_gene_tree_newick\n";
  #   $first = 0;
  # }
  my $orthologger_outstring = $orthologger_obj->ortholog_result_string();

  store_orthologger_out($orthologger_outstring, [\%idpair_orthocount_all, \%idpair_orthocountMB], $MB_count);
  $orthologger_obj->decircularize();
}
}

# output

open my $fh_ortho, ">$output_orthologs_filename";

$n_bootstrap = max($n_NJ_bootstrap, $n_ML_bootstrap);
my $iwidth = length(max($n_rep, $n_bootstrap));
my $iformat = '%' . $iwidth . 'i';
my $id1_old = '';
my @sorted_idpairs = sort {$a cmp $b} keys %idpair_orthocount_all;
foreach my $idpair (@sorted_idpairs) {
  my $actual_countNJ = (defined $idpair_actual_orthocountNJ{$idpair})? 
    $idpair_actual_orthocountNJ{$idpair} : 0;
  my $actual_countML = (defined $idpair_actual_orthocountML{$idpair})? 
    $idpair_actual_orthocountML{$idpair} : 0;
  my $bs_countNJ = (defined $idpair_bs_orthocountNJ{$idpair})? $idpair_bs_orthocountNJ{$idpair} : 0;
  my $n_bootstrap_ML = ($tree_find_method eq 'ML')? $n_ML_bootstrap : 0;
  my $n_rep_ML = ($tree_find_method eq 'ML')? $n_rep : 0;
  my $bs_countML = (defined $idpair_bs_orthocountML{$idpair})?  $idpair_bs_orthocountML{$idpair} : 0;
  my $count_MB = (defined $idpair_orthocountMB{$idpair})? $idpair_orthocountMB{$idpair} : 0;
  next if( ($actual_countNJ == 0) and
	   ($actual_countML == 0) and
	   ($bs_countNJ/$n_NJ_bootstrap) < $min_bs_support and
	   ($bs_countML/$n_ML_bootstrap) < $min_bs_support and
	      ($count_MB/$sum_MB_count) < $min_bs_support);
  my ($id1, $id2) = split(" ", $idpair);
  if ($id1 ne $id1_old) {
    printf($fh_ortho "\n%-38s NJ($iformat)  ML($iformat) bootstrap: NJ($iformat)  ML($iformat)  MB($iformat)\n", 
	   $id1, $n_rep, $n_rep_ML, $n_NJ_bootstrap, $n_bootstrap_ML, $sum_MB_count);
    $id1_old = $id1;
  }
  printf($fh_ortho "    %-37s $iformat      $iformat                $iformat      $iformat    $iformat\n", 
	 $id2, $actual_countNJ, $actual_countML, $bs_countNJ, $bs_countML, $count_MB);
}
close $fh_ortho;

sub run_clearcut{
  my $overlap_fasta_string = shift;
  my $clearcut_cl = shift;
  open my $fhtmp, ">tmp_alignment_fasta";
  print $fhtmp $overlap_fasta_string, "\n";
  close $fhtmp;
  #  print STDERR "about to run clearcut. $clearcut_cl \n";
  my $clearcut_newicks_out = `$clearcut_cl < tmp_alignment_fasta`;
  # print STDERR "clearcut finished.\n";
  return $clearcut_newicks_out;
}

sub run_quicktree{
  my $overlap_fasta_string = shift;
  my $correction = shift || 'kimura';
  open my $fhtmp, ">tmp_overlap_fasta";
  print $fhtmp $overlap_fasta_string, "\n";
  close $fhtmp;

  my $overlap_stockholm_string = `sreformat stockholm tmp_overlap_fasta > tmp_overlap_stockholm`;
  my $newick_out;
  if ($correction eq 'kimura') {
    $newick_out = `quicktree -kimura tmp_overlap_stockholm`;
  } else {
    $newick_out = `quicktree  tmp_overlap_stockholm`;
  }

  $newick_out =~ s/\s+//g;
  return $newick_out;
}

sub run_fasttree{
  my $overlap_fasta_string = shift;
  my $fasttree_cl = shift;
  open my $fhtmp, ">tmp_alignment_fasta";
  print $fhtmp $overlap_fasta_string, "\n";
  close $fhtmp;
  #  print STDERR "about to run fasttree. $fasttree_cl \n";
  my $fasttree_newick_out = `$fasttree_cl  tmp_alignment_fasta`;
  # print STDERR "fasttree finished.\n";
  return $fasttree_newick_out;
}


sub store_orthologger_out{
  my $orthologger_outstring = shift;
  my $idpair_count_hrefs = shift;
  my $weight = shift || 1;
  my @orthologger_out = split("\n", $orthologger_outstring);
  foreach my $id_orthologs (@orthologger_out) {
    if ($id_orthologs =~ /orthologs of\s+(\S+):\s+(.*)/) {
      my $id1 = $1;
      my @orthologs = split(" ", $2);
      foreach my $id2 (@orthologs) {
	my $idpair = $id1 . "   " . $id2;
	foreach my $idpair_orthocount (@$idpair_count_hrefs) {
	  $idpair_orthocount->{$idpair} += $weight;
	}
      }	      # loop over ortholog names in line of orthologger output
    }
  }			       # loop over lines in orthologger output
}


# Bayesian
# we have set of gene tree newicks (for topologies visited by markov chain)
# each one has associated vector of hits in the various runs done.
# Want to analyze each topology with Orthologger, for each potential
# ortholog relationship, keep track of the number of hits (just sum over runs for now).

sub process_MB_out{
  my $alignment_nex_filename = shift; # e.g. fam9877.nex
  my $burnin_frac = shift || 0.1;

  my $mrb_obj = Mrbayes->new({'alignment_nex_filename' =>$alignment_nex_filename});
  my $newick_count = {}; # ref to hash of (newick(with ids): hit_count) pairs
	# check if MrBayes output files are present:
	my $MB_output_present = 1;
	my $run_filename = $alignment_nex_filename . "run1.t";
	if(! (-f $alignment_nex_filename and -f $run_filename)){
	warn "MrBayes output not available, no Bayesian analysis performed.\n";
	return $newick_count;
	}

  my $gen_param_hrefs = 
    #load_params($alignment_nex_filename);
    $mrb_obj->retrieve_param_samples(); # read in from *.run?.p file
  my ($gen_ntopo_hrefs, $newick_number_map, $number_newick_map) = 
    #  load_topologies($alignment_nex_filename);	#read in from * .run?.t file
    $mrb_obj->retrieve_topology_samples();
  # $topo_count is
  my ($topo_count, $total_trees) = 
    #  count_topologies($gen_ntopo_hrefs);
    $mrb_obj->count_topologies($gen_ntopo_hrefs);
  my $distinct_newicks = scalar keys %$topo_count;
  print "Distinct topologies: $distinct_newicks\n";

  my @sorted_tree_numbers = sort { sum(@{$topo_count->{$b}}) <=> sum(@{$topo_count->{$a}}) } keys %$topo_count; # tree numbers, sorted by occurrences (sum of runs).

  my $number_id_map = 
    #  load_number_id_map($alignment_nex_filename);
    $mrb_obj->retrieve_number_id_map();
  print "n trees (sorted): ", scalar @sorted_tree_numbers, "\n";
  my  %number_rank_map = (); # number is n for nth distinct topology visited by chain
  my %rank_number_map = (); # rank is 1 for topology with highest posterior prob. etc.
  my $index = 1;
  my $total_trees_so_far = 0;
  my $sum_diff = 0;
  my $L1_denom = 0; 
  foreach my $tree_number (@sorted_tree_numbers) {
    my $treecount = sum @{$topo_count->{$tree_number}};
    my @treecounts = @{$topo_count->{$tree_number}};
    $L1_denom += $treecount;
    $sum_diff += max(@treecounts) - min(@treecounts);
    if ($treecount % 2 == 1) {
      $sum_diff--;
      $L1_denom--;
    }

    #  print "treecounts: ", join(" ", @treecounts), "\n";
    $total_trees_so_far += $treecount;
    my $newick_with_numbers = $number_newick_map->{$tree_number};
    my $post_prob = $treecount/$total_trees;
    $number_rank_map{$tree_number} = $index;
    $rank_number_map{$index} = $tree_number;
    printf("%4i %4i %4i %4i %4i  %10.5g %10.5g  %4i  %s\n", 
	   $index, @treecounts, $treecount, $total_trees_so_far, 
	   $post_prob, $total_trees_so_far/$total_trees, $tree_number,
	   $newick_with_numbers);
    my $newick_with_ids = $mrb_obj->restore_ids_to_newick($newick_with_numbers, $number_id_map);

    $newick_count->{$newick_with_ids} = $treecount;
    $index++;
  }
  foreach (keys %$newick_count){
print "count, newick: ", $newick_count->{$_}, " $_ \n";
}
  #print "$sum_diff, $total_trees $L1_denom\n";
  printf("topology post. distrib. L1 difference between runs: %8.5f \n", 0.5*$sum_diff/$L1_denom);
  print "total tree hits: $total_trees \n";
  print "distinct newicks: $distinct_newicks \n";

  my $j_run = 1;
  foreach my $gen_ntopo (@$gen_ntopo_hrefs) {
    my $string = "Run: $j_run\n";
    my $trees_read_in = scalar keys %{$gen_ntopo};
    print "run, trees read in $j_run  $trees_read_in \n";
    my $n_burnin = int($burnin_frac * $trees_read_in);
    my @sorted_generations = sort {$a <=> $b} keys %{$gen_ntopo};
    # print "nburnin: $n_burnin  ", scalar @sorted_generations, "\n";
    foreach my $i_gen (@sorted_generations[$n_burnin..$trees_read_in-1]) {
      $string .= "$i_gen " .
	$number_rank_map{$gen_ntopo->{$i_gen}} . " " .
	  $gen_param_hrefs->[$j_run-1]->{$i_gen} . "\n";
    }				# loop over gens
    #  print "$string\n" if(1);
    open my $fhtp, ">run$j_run.tp";
    print $fhtp "$string\n";
    close $fhtp;
    $j_run++;
  }
  return $newick_count;
}

# (jgi_Selmo1_134275,jgi_Selmo1_448960,(((((((evm.model.supercontig_3.375,(Glyma05g33720.1,Glyma08g06000.1)),(Bradi4g29810.1,(((Sb07g009430.1,POPTR_0012s04240.1),POPTR_0015s04860.1),X_29653.m000293))),GSVIVT01008456001),GRMZM2G357034_P01),LOC_Os09g23640.1),Solyc01g097430.2.1),IMGA_Medtr8g107450.1))
