#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw ( min max sum );

use lib '/home/tomfy/Orthologger/lib';
use Overlap;
use Orthologger;

use lib '/home/tomfy/cxgn/cxgn-corelibs/lib';
use CXGN::Phylo::File;
use CXGN::Phylo::Parser;

use Devel::Cycle; # for finding circular refs - cause of memory leaks.
# find_cycle($test); # to find circular refs in $test (which might be an object, e.g.)

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

# print "# $0 ", join(" ", @ARGV), "\n";

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
# -m <min bootstrap support> (Output only orthologs having >= this bootstrap support; if > 1 interpret as %)

# typical usage: perl bootstrap_ortholog.pl -i $align_filename -T ML -k -n 1 -N 100 -S 12345 -r mindl -m 0.015 -q castorbean > outfile

# get options
getopts("i:s:S:f:T:N:q:Q:kn:r:m:");

die "usage: bootstrap_ortholog.pl -i <input_filename> [-s species_tree -S <rng seed> -f <overlap fraction> -T <tree infer type>  -N <#bootstrap data sets> -q <query sequence id> -k -n <# rng trees from each data set> i" unless($opt_i and -f $opt_i);

my $input_file = $opt_i;

my $n_bootstrap = $opt_N || 0; # default is just do the actual data, not bootstraps.
my $nongap_fraction = $opt_f || 0.8;

my $query_species = ($opt_q)? $opt_q: undef;
my $query_id_regex = ($opt_Q)? $opt_Q: undef;

my $bootstrap_seed = ($opt_S)? $opt_S: undef;
my $tree_find_seed = ($opt_S)? $opt_S + $seed_increment: undef;

my $clearcut_base_cl = 'clearcut --stdin --stdout -a -P ';
#$clearcut_base_cl = 'MrRogers --stdin --stdout -a -P ';
my $fasttree_base_cl = 'FastTree -wag -gamma -bionj -quiet ';

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
my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this.
$align_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg; # to make clearcut happy.
#exit;
# construct an overlap object.
my $overlap_obj = Overlap->new($align_string, $nongap_fraction, $bootstrap_seed);
my $mb_output_string = '';
my $seed = 1234567;
my $swapseed = 135786421;
my $n_temperatures = 4;
my $n_runs = 2;
my $temperature_gap = 0.25;	# 333; # temperature spacing of chains
my $chunk_size = 100;	   # steps between decisions whether to go on.
my $print_freq = int($chunk_size/2);
my $sample_freq = 10;		# int($chunk_size/100);

$sample_freq = max($sample_freq, 1);
$print_freq = max($print_freq, 1);
my $mrb_outfile_basename = $input_file;

$mrb_outfile_basename =~ s/[.]fasta//;
$mrb_outfile_basename .= '.nex';
my $mc3swapstats_filename = $mrb_outfile_basename . ".mc3stat";
open my $fhmc3, ">$mc3swapstats_filename";

open my $fh1, ">$mrb_outfile_basename";
my $overlap_nexus_string = $overlap_obj->overlap_nexus_string();
print $fh1 $overlap_nexus_string, "\n";
close $fh1;

my $ngen= $chunk_size;
my $n_burnin_samples = int(0.25*$ngen/$sample_freq);
my $burnin_frac = 0.125;
my $diag_freq = $chunk_size;	# int($chunk_size/2);
# Allchains       Yes/No                No
#  Allcomps        Yes/No                No
#  Relburnin       Yes/No                Yes
#  Burnin          <number>              0
#  Burninfrac      <number>              0.25
my $begin_piece =  "begin mrbayes;\n" .
  "set autoclose=yes nowarn=yes;\n";
my $seed_piece = '';
if(defined $seed){ $seed_piece .=  "set seed=$seed;\n"; } 
if(defined $swapseed){ $seed_piece .= "set swapseed=$swapseed;\n"; }
my $middle_piece =  "execute $mrb_outfile_basename;\n" .
  #  "set scientific=no;\n" .
  "set precision=6;\n" .
  "lset rates=invgamma;\n" .
  "prset aamodelpr=fixed(wag);\n" .
  # "prset pinvarpr=fixed(0.15);\n" .
  #"mcmcp mcmcdiagn=no;\n" .
  #"mcmcp diagnfreq=$diag_freq;\n" .
  #"mcmcp diagnstat=maxstddev;\n" .
  "mcmcp minpartfreq=0.02;\n" . # bipartitions with freq. less than this are not used in the  diagnostics (default is 0.10)
  #"mcmcp allchains=yes;\n" .
  #"mcmcp allcomps=yes;\n" .
  #"mcmcp relburnin=no;\n" .
  "mcmcp burninfrac=$burnin_frac;\n" .
  #  "mcmcp ordertaxa=yes;\n" .
  "mcmcp nchains=$n_temperatures;\n" .
  "mcmcp nruns=$n_runs;\n" .
  "mcmcp temp=$temperature_gap;\n" .
  "mcmcp samplefreq=$sample_freq;\n" .
  "mcmcp printfreq=$print_freq;\n" . 
  "mcmcp checkpoint=yes checkfreq=$chunk_size;\n";

my $summary_piece = "sump;\n" . "sumt;\n";
my $end_piece = "end;\n";

my $mrbayes_nexus1 = $begin_piece .
  $seed_piece .
  $middle_piece .
  "mcmc ngen=$ngen;\n" .
  $summary_piece . 
  $end_piece;

open my $fh, ">tmp_mrb1.nex";
print $fh "$mrbayes_nexus1";
close $fh;
#print "about to run mb in backticks\n";
$mb_output_string .= `mb tmp_mrb1.nex`;
#print "after mb in backticks\n";
#print $mb_output_string, "\n";
#system "cat $mrb_outfile_basename.tstat";



# my ($n_splits_used, $avg_stddev, $splits_n_bad) = splits_convergence($mrb_outfile_basename);
# my ($mpc_string, $mp_n_bad) = modelparam_convergence($mrb_outfile_basename);
# print "$n_splits_used $avg_stddev $splits_n_bad   ";
# print "$mpc_string $mp_n_bad \n";
#exit;

print $fhmc3 "$ngen ", extract_swap_info($mb_output_string), "\n";

open $fh, ">mb.stdout";
print $fh "$mb_output_string \n";
close $fh;

my ($converge_count, $conv_string) = test_convergence($mrb_outfile_basename);
print "0 $ngen $converge_count  $conv_string\n";

foreach my $i (1..100000) {
  $ngen += $chunk_size;
  #$n_burnin_samples = int(0.25*$ngen/$sample_freq);
  my $mrbayes_nexus2 = $begin_piece .
    $middle_piece .
      "mcmc ngen=$ngen burninfrac=$burnin_frac append=yes;\n" .
	"sump;\n" . "sumt;\n" .
	  $end_piece;

  open $fh, ">tmp_mrb2.nex";
  print $fh "$mrbayes_nexus2";
  close $fh;

  $mb_output_string =  `mb tmp_mrb2.nex`;

  print $fhmc3 "$ngen ", extract_swap_info($mb_output_string), "\n";
  open $fh, ">mb2.stdout";
  print $fh "$mb_output_string \n";
  close $fh;

  my $conv_count;
  ($conv_count, $conv_string) = test_convergence($mrb_outfile_basename);
  $converge_count += $conv_count;
  print "$i $ngen $converge_count  $conv_string\n";
  last if($converge_count > 20);
}

exit;

$ngen += $chunk_size;
# my $mrbayes_nexus3 =
# $begin_piece .
#   $middle_piece .
#   "mcmc ngen=$ngen append=yes;\n" .
#   $summary_piece .
#   $end_piece;
my $mrbayes_nexus3 = "begin mrbayes;\n" .
  "set autoclose=yes nowarn=yes;\n" .
  "execute $mrb_outfile_basename;\n" .
  "mcmcp burninfrac=$burnin_frac;\n" .
  "sump;\n" .
  "sumt;\n" .
  "end;\n";

open $fh, ">tmp_mrb3.nex";
print $fh "$mrbayes_nexus3";
close $fh;

$mb_output_string = `mb tmp_mrb3.nex`;
open $fh, ">mb3.stdout";
print $fh "$mb_output_string \n";
close $fh;
exit;


sub splits_convergence{
  my $file_basename = shift;	# e.g. fam9877.nex
  my $max_ok_stddev = shift || 0.1; # convergence is 'ok' for a split if stddev < this.
  my $min_hits = shift || 25; # ignore splits with fewer hits than this.

  my $filename = $file_basename . ".tstat";
  open my $fh, "<$filename";
  my @lines = <$fh>;

  my ($avg_stddev, $count, $bad_count) = (0, 0, 0);
  foreach (@lines) {
    #   print;
    next unless(/^\s*\d/); # skip if first non-whitespace is not numeral.
    my @cols = split(" ", $_);
    my $hits = $cols[1]; my $stddev = $cols[3];
    #  print "$hits, $min_hits, $stddev\n";
    last if($hits < $min_hits);
    $count++;
    $avg_stddev += $stddev;
    if ($stddev > $max_ok_stddev) {
      $bad_count++;
      next;
    }
  }
  $avg_stddev = ($count == 0)? 100 : $avg_stddev/$count;
  return ($count, $avg_stddev, $bad_count);
}


sub modelparam_convergence{
  my $file_basename = shift;
  my $min_ok_ESS = shift || 100;
  my $max_ok_PSRF = shift || 1.1;
  my $n_bad = 0;
  my $string = '';

  open $fh, "<$file_basename.pstat";
  my @lines = <$fh>;
  my $discard = shift @lines;
  foreach (@lines) {
    my @cols = split(" ", $_);
    my ($minESS, $PSRF) = @cols[6,8];
    next unless($minESS =~ /^\d*[.]?\d+/ and $PSRF =~ /^\d*[.]?\d+/);
    $string .= "$minESS $PSRF ";
    if ($minESS < $min_ok_ESS  or  $PSRF > $max_ok_PSRF) {
      $n_bad++;
    }

  }
  return ($string, $n_bad);
}

sub test_convergence{
  my $file_basename = shift;	# e.g. fam9877.nex

  my $splits_max_ok_stddev = shift || 0.1; # convergence is 'ok' for a split if stddev < this.
  my $splits_min_hits = shift || 25; # ignore splits with fewer hits than this.

  my $model_param_min_ok_ESS = shift || 100;
  my $model_param_max_ok_PSRF = shift || 1.1;

  my ($splits_count, $splits_avg_stddev, $splits_bad_count) = 
    splits_convergence($file_basename, 
		       $splits_max_ok_stddev, $splits_min_hits);
  my ($modelparam_string, $modelparam_n_bad) = 
    modelparam_convergence($file_basename, 
			   $model_param_min_ok_ESS, $model_param_max_ok_PSRF);

  my $conv_string = "$splits_count $splits_avg_stddev $splits_bad_count " . 
    " $modelparam_string  $modelparam_n_bad ";
  my ($ESS1, $PSRF1, $ESS2, $PSRF2) = split(" ", $modelparam_string);
  #  my $bow1 = ($PSRF1-1); # bow = psrf**2-1
  #  my $bow2 = ($PSRF2-1);
  my $converged = (($splits_bad_count < 1) and ($splits_avg_stddev < 0.01) and
		   ($ESS1 > 100 and $ESS2 > 100) and ($PSRF1 < 1.05 and $PSRF2 < 1.05));
  return ($converged? 1 : 0, $conv_string);
}


sub extract_swap_info{
  my $mb_stdout_string = shift;
  my @mb_stdout_lines = split("\n", $mb_stdout_string);
  my $n_lines_to_extract = 0;
  my $extract_next_n = 0;
  my $n_runs = undef;
  my $n_chains = undef;
  my $out_string = '';
  foreach (@mb_stdout_lines) {
    if (/number of chains to (\d+)/) {
      $n_chains = $1;
      $n_lines_to_extract = $n_chains + 4;
      last if(defined $n_runs);
    } elsif (/number of runs to (\d+)/) {
      $n_runs = $1;
      last if(defined $n_chains);
    }

  }
  my $run;
  my %run_string = ();
  foreach (@mb_stdout_lines) {
    if (/Chain swap information for run (\d+)/) {
      $run = $1;
      $extract_next_n = $n_lines_to_extract;
    }
    $out_string .= "$_\n" if($extract_next_n > 0);
    $extract_next_n--;
    if ($extract_next_n == 0) {
      $run_string{$run} = $out_string;
      $out_string = '';
      last if($run == $n_runs);
    }
  }
  $out_string = '';

  foreach (keys %run_string) {
    #print "$_  ", $run_string{$_}, "\n";
    my @lines = split("\n", $run_string{$_});
    splice @lines, 0, 4;
    #  print join("\n", @lines);

    my %ij_swap_pA = ();
    my %ij_swap_tries = ();
    foreach my $i (1..$n_chains) {
      my $l = $lines[$i-1];
      $l =~ s/^\s*\d+\s+[|]\s+//;
      my @xs = split(" ", $l);
      my $n_ntry = $i-1;
      my $n_pA = $n_chains-$i;
      foreach my $j (1..$n_ntry) {
#	print "swap_tries key: [$i $j]\n";
	$ij_swap_tries{"$i $j"} = shift @xs;
      }
      foreach (1..$n_pA) {
	my $j = $_ + $i;
#	print "swap_pA key: [$j $i]\n";
	$ij_swap_pA{"$j $i"} = shift @xs;
      }
    }				# loop over chains
    my %ij_swap_accepts = ();
    # my @sijs = sort {$a cmp $b} keys %ij_swap_tries;
    # foreach (@sijs) {

    foreach my $diff (1..$n_chains-1) {
      foreach my $i (1..$n_chains-1) {
	my $j = $i + $diff;

	last if($j > $n_chains);
	my $key = "$j $i";
#	print "i,j: $i, $j, key: [$key] \n";
	if (exists $ij_swap_pA{$key} and exists $ij_swap_tries{$key}) {
	  $ij_swap_accepts{$key} = $ij_swap_tries{$key} * $ij_swap_pA{$key};
	  $out_string .= int($ij_swap_accepts{$key}+0.5) . " " . $ij_swap_tries{$key} . "  ";
	} else {
	  warn "key $key present in neither ij_swap_tries nor ij_swap_pA.\n";
	}
      }
      $out_string .= ' ';
    }
    $out_string .= ' ';
  }				# loop over runs
  return $out_string;
}
