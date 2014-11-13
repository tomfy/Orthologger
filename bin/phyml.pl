#!/usr/bin/perl -w
use strict;

use lib '/home/tomfy/Orthologger/lib';
use Phyml;

use Getopt::Std;

use vars qw($opt_i $opt_o $opt_b $opt_u $opt_v $opt_a $opt_m $opt_c $opt_d);
# -i <fasta_alignment_file>
# -d <'dna' or 'protein'> 
# -o <what to optimize, l, t, r, lr, tl, tlr>
# -b <bootstraps. n>1 -> number of bootstrap replicates. n<0 various support tests. n=0 nothing>
# -u <initial tree newick file>
# -v <proportion invariant sites. e.g.: -v 0.08  or  -v e  to optimize pinv>
# -a <discrete gamma model shape parameter alpha. e.g.: -a 0.6   or  -a e to optimize alpha>
# -s <NNI or SPR or BEST> # topological moves to try.

# get options
getopts("i:o:b:u:v:a:m:c:d:");

my $dna_or_protein = (defined $opt_d)? $opt_d : 'protein';
my $align_fasta_file = (defined $opt_i)? $opt_i : die "no input alignment file specified.\n";
my $optimize_param = (defined $opt_o)? $opt_o : 'none';
my $n_bootstrap = (defined $opt_b)? $opt_b : 0;
my $initial_tree_newick_file = (defined $opt_u)? $opt_u : undef;
# opt_v, $opt_a not implemented yet.
my $alpha = (defined $opt_a)? $opt_a : undef;
my $n_rate_classes = (defined $opt_c)? $opt_c : undef;
my $p_invariant = (defined $opt_v)? $opt_v : undef;
my $subst_model = (defined $opt_m)? $opt_m : undef;

# print "p_invariant: [$p_invariant] \n";

my $phyml_obj = Phyml->new({
			    'dna_or_protein' => $dna_or_protein,
			    'fasta_file' => $align_fasta_file,
			   'optimize_param' => $optimize_param,
			    'n_bootstrap' => $n_bootstrap,
			    'initial_tree_newick_file' => $initial_tree_newick_file,
			    'alpha' => $alpha,
			    'p_invariant' => $p_invariant,
			    'subst_model' => $subst_model,
			    'n_rate_classes' => $n_rate_classes
});

$phyml_obj->run();

print "# phyml command line: ", $phyml_obj->{phyml_command_line}, "\n";
print $phyml_obj->{phyml_stdout}, "\n";

print "# ln(L): ", $phyml_obj->{log_likelihood}, "\n";
# print "# tree: ", $phyml_obj->{newick_out}, "\n";
