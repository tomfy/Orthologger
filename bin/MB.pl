#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Std;
use Getopt::Long;
use List::Util qw ( min max sum );

#use lib '/home/tomfy/Orthologger/lib';
use CXGN::Phylo::Overlap;
use CXGN::Phylo::Mrbayes; # perl module encapsulating mrbayes bayesian phylogeny program.

#use Devel::Cycle; # for finding circular refs - cause of memory leaks.
# find_cycle($test); # to find circular refs in $test (which might be an object, e.g.)

# read in an alignment file. get an overlap object. 

#use vars qw($opt_i $opt_S $opt_f);
# -i <input_alignment_filename>     (name of input alignment file, fasta format)
# -S <seed>  (rng seed, default: get from clock)
# -f <fraction> (fraction of non-gap chars required in alignment column to keep it. Default is 0.8)


# typical usage: perl MB.pl -i fam9877.nex

# get options
#getopt("i:S:f:");

# defaults:
my $input_file = undef;
my $seed = undef;
my $nongap_fraction = 0.8;
my $chunk_size = 200;
my $n_temperatures = 4;
my $delta_temperature = 0.1;
my $sample_freq = 20;
my $n_runs = 3;
my $burnin_fraction = 0.1;
my $converged_chunks_required = 10;
my $modelparam_min_ok_ESS = 200;


GetOptions('input_file=s' => \$input_file,
	   'seed=i' => \$seed,
	   'nongap_fraction=f' => \$nongap_fraction,
	   'chunk_size=i' => \$chunk_size,
	   'n_temperatures=i' => \$n_temperatures,
	   'delta_temperature=f' => \$delta_temperature,
	   'sample_freq=i' => \$sample_freq,
	   'n_runs=i' => \$n_runs,
	   'burn-in_fraction=f' => \$burnin_fraction,
	   'converged_chunks_required=i' => \$converged_chunks_required,
	   'ESS_min=i' => \$modelparam_min_ok_ESS);
print "Seed: $seed\n";
print "chunksize: $chunk_size \n";
print "nongapfrac: $nongap_fraction\n";
print "n_temperatures: $n_temperatures\n";
print "delta temperature: $delta_temperature \n";
print "sample_freq: $sample_freq\n";
print "n_runs: $n_runs\n";
print "burnin fraction: $burnin_fraction\n";
print "converged chunks required: $converged_chunks_required\n";
print "modelparam min OK ESS: $modelparam_min_ok_ESS\n";

# exit;

die "Must specify name of input alignment file. Input file: [$input_file].\n" . "Usage: MB.pl --input_file <input_filename> [--seed <rng seed> --nongap_fraction <overlap fraction>]\n" unless($input_file and -f $input_file);

#my $input_file = $opt_i;


#### Get the alignment:
my $align_string =  `cat $input_file`;
# fixes to $align_string:
$align_string =~ s/IMGA[|]/IMGA_/g; #pipes in id cause problem; replace '|' with '_'.
my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this.
$align_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg; # to make clearcut happy.

# construct an overlap object.
# my $bootstrap_seed = 1234567;	# ($opt_S)? $opt_S : undef;
#my $nongap_fraction = ($opt_f)? $opt_f : 0.8;
my $overlap_obj = CXGN::Phylo::Overlap->new($align_string, $nongap_fraction); # , $bootstrap_seed);

# construct MrBayes object and run
#my $seed = ($opt_S)? $opt_S : undef;
my $swapseed = ($seed)? $seed + 1000 : undef;
#my $n_temperatures = 4;
my $n_swaps = ($n_temperatures > 2)? int ( ($n_temperatures-1) * $n_temperatures/2 / 2.1 ) : 1; # 1,2,3 T's -> 1, 4 T's ->2, 5 T's -> 4
print "n_swaps: $n_swaps \n";
#my $n_runs = 2;
#my $temperature_gap = 0.3;   # temperature spacing of chains
#my $chunk_size = 200;	   # steps between decisions whether to go on.
#my $print_freq = int($chunk_size/2);
#my $sample_freq = 20;
$sample_freq = max($sample_freq, 1);


my $mrb_outfile_basename = $input_file;

$mrb_outfile_basename =~ s/[.]fasta//;
my $alignment_nex_filename = $mrb_outfile_basename . '.nex';

open my $fh1, ">$alignment_nex_filename";
my $overlap_nexus_string = $overlap_obj->overlap_nexus_string();

print $overlap_nexus_string, "\n";
print $fh1 $overlap_nexus_string, "\n";
close $fh1;

my $mrb_obj = CXGN::Phylo::Mrbayes->new(
			   {'alignment_nex_filename' =>$alignment_nex_filename,
			    'chunk_size' => $chunk_size,
			    'seed' => $seed,
			    'swapseed' => $swapseed,
			    #	   'fixed_pinvar' => 0.4
			    'sample_freq' => $sample_freq,
			    'modelparam_min_ok_ESS' => $modelparam_min_ok_ESS,
			    'n_temperatures' => $n_temperatures,
			    'delta_temperature' => $delta_temperature,
			    'n_swaps' => $n_swaps,
			    'n_runs' => $n_runs,
			    'burnin_fraction' => $burnin_fraction,
			    'converged_chunks_required' => $converged_chunks_required,
#			    'temperature_factor' => 1.414, # I wanted to have T's exponentially spaced, but mb does not allow
			   }
			  );
$mrb_obj->run();

exit;
