#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw ( min max sum );

#use lib '/home/tomfy/Orthologger/lib';
use CXGN::Phylo::Overlap;
use CXGN::Phylo::Mrbayes; # perl module encapsulating mrbayes bayesian phylogeny program.

#use Devel::Cycle; # for finding circular refs - cause of memory leaks.
# find_cycle($test); # to find circular refs in $test (which might be an object, e.g.)

# read in an alignment file. get an overlap object. 

use vars qw($opt_i $opt_S $opt_f);
# -i <input_alignment_filename>     (name of input alignment file, fasta format)
# -S <seed>  (rng seed, default: get from clock)
# -f <fraction> (fraction of non-gap chars required in alignment column to keep it. Default is 0.8)


# typical usage: perl MB.pl -i fam9877.nex

# get options
getopts("i:S:f:");

die "usage: MB.pl -i <input_filename> [-S <rng seed> -f <overlap fraction>\n" unless($opt_i and -f $opt_i);

my $input_file = $opt_i;


#### Get the alignment:
my $align_string =  `cat $input_file`;
# fixes to $align_string:
$align_string =~ s/IMGA[|]/IMGA_/g; #pipes in id cause problem; replace '|' with '_'.
my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this.
$align_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg; # to make clearcut happy.

# construct an overlap object.
# my $bootstrap_seed = 1234567;	# ($opt_S)? $opt_S : undef;
my $nongap_fraction = ($opt_f)? $opt_f : 0.8;
my $overlap_obj = Overlap->new($align_string, $nongap_fraction); # , $bootstrap_seed);

# construct MrBayes object and run
my $seed = ($opt_S)? $opt_S : undef;
my $swapseed = ($opt_S)? $opt_S + 1000 : undef;
#my $n_temperatures = 4;
#my $n_runs = 2;
#my $temperature_gap = 0.3;   # temperature spacing of chains
my $chunk_size = 1000;	   # steps between decisions whether to go on.
#my $print_freq = int($chunk_size/2);
my $sample_freq = 20;
$sample_freq = max($sample_freq, 1);

my $modelparam_min_ok_ESS = 500;

my $mrb_outfile_basename = $input_file;

$mrb_outfile_basename =~ s/[.]fasta//;
my $alignment_nex_filename = $mrb_outfile_basename . '.nex';

open my $fh1, ">$alignment_nex_filename";
my $overlap_nexus_string = $overlap_obj->overlap_nexus_string();

print $overlap_nexus_string, "\n";
print $fh1 $overlap_nexus_string, "\n";
close $fh1;

my $mrb_obj = Mrbayes->new(
			   {'alignment_nex_filename' =>$alignment_nex_filename,
			    'chunk_size' => $chunk_size,
			    'seed' => $seed,
			    'swapseed' => $swapseed,
			    #	   'fixed_pinvar' => 0.4
			    'sample_freq' => $sample_freq,
			    'modelparam_min_ok_ESS' => $modelparam_min_ok_ESS,
			   }
			  );
$mrb_obj->run();

exit;
