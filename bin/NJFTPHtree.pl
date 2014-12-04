#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
use List::Util qw (min max sum);


no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

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

use Phyml;
use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info';

use CXGN::Phylo::Overlap;
use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

# The idea here is to get several NJ bootstrap trees, and use them as initial trees
# for FastTree to optimize (w.r.t. actual data, not the bootstrap resampled data).
# and then get the few with the best likelihoods, and further optimize using Phyml.

my $gg_filename;  # required  cl parameter 
# by default do just fasttree, and no bootstrapping.
my $n_bs = 5;
#my ($do_ft, $n_ft_bs) = (1, 0);
my ($do_phyml, $n_phyml_bs) = (0, 0);
my $phyml_opt = 'r';

my $input_alignment_file = undef;
my $output_newick_file = undef;
# print "[$gg_filename] [$do_nj] [$do_ml] [$n_bs] [$ml_bs] \n";
# my $n_taxa                   = shift || 21;
my $species_tree_newick_file = undef;
my $reroot_method            = 'mindl';
my $nongap_fraction = 0.8;
my $min_overlap_length = 40;
my $max_overlap_length = 1000000000;
my $min_n_sequences = 0;
my $max_n_sequences = 1000000;
my $support = 1;		# default is to do support
my $additional_ft_options = '';



# At present:
# protein alignment only as input
# WAG
# for phyml, 4 rate categories plus invariant, alpha and p_inv estimated
GetOptions(
	   'gg_file=s'           => \$gg_filename, #
	   #	   'nj!'          => \$do_nj, # whether to do NJ tree construction for actual data
	   #  exclamation point means can use  -nj  and  -nonj
	   #	   'ft!'      => \$do_ft, # whether to do FastTree for actual data
	   'phyml!' => \$do_phyml, # whether to do Phyml tree construction for actual data
	   'phyml_opt=s' => \$phyml_opt,
	   #	   'support!'          => \$support, # -nosupport to turn off outputting of local branch support numbers.
	   'n_bs=i' => \$n_bs,	# number of NJ bootstraps to do
	   #	   'n_ft_bs=i' => \$n_ft_bs, # number of FastTree bootstraps to do
	   #	   'n_phyml_bs=i' => \$n_phyml_bs, # number of Phyml bootstraps to do
	   #    'ml_bs!' => \$ml_bs, # boolean to do ML bs or not (default: 0)
	   #	   'species_tree=s' => \$species_tree_newick_file,
	   #	   'reroot_method=s' => \$reroot_method, 
	   'nongap_fraction=f' => \$nongap_fraction,
	   'min_overlap_length=i' => \$min_overlap_length,
	   'max_overlap_length=i' => \$max_overlap_length,
	   'max_n_sequences=i' => \$max_n_sequences,
	   'input_alignment_file=s' => \$input_alignment_file,
	   'output_newick_file=s' => \$output_newick_file,
	   #	   'ft_options=s' => \$additional_ft_options,
	  );
if (!defined $gg_filename  or  ! -f $gg_filename) {
  # die "No gene-genome association file specified. Exiting. \n";
}
my $fh_in = \*STDIN;
if (defined $input_alignment_file) {
  open $fh_in, "<", "$input_alignment_file" or die "File [$input_alignment_file] could not be opened for reading.\n";
}
my $fh_out;
if (defined $output_newick_file) {
  open $fh_out, ">", "$output_newick_file" or die "File [$output_newick_file] could not be opened for writing.\n";
}else{
 $fh_out = \*STDOUT;
}

print STDERR # "gg_file: ", (defined $gg_filename)? $gg_filename: 'undefined', " \n", 
  "n_bs: $n_bs \n", 
  "do_phyml: $do_phyml \n",
  "min nongap_fraction: $nongap_fraction \n",
  "min_overlap_length: $min_overlap_length \n";
my $state              = 'idline'; # other possible values: fasta
my ( $qid, $fam_size, $taxa, $idline, $fasta);
my $support_string = ($support)? '' : '-nosupport'; # controls whether FastTree calculates local branch support numbers.

my ($n_ft, $n_phyml) = (0, 0);

while (<$fh_in>) {
  if ( /^Id\s+(\S+)\s+family/ ) {
    $qid = $1;
    my @cols = split( " ", $_ );
    ( $qid, $fam_size, $taxa) = @cols[ 1, 4, 5];
    $idline = $_;
    $state  = 'fasta';
    $fasta  = '';
  } elsif ( $state eq 'fasta' ) {

    if (/^\s*$/) { # blank line after sequence -> process the sequence.
      chomp $idline;
      my $string_to_print = "$idline   ";
      my $rng_seed = srand();
      if ($fasta ne '') {
	my $overlap_obj =
	  CXGN::Phylo::Overlap->new( $fasta, $nongap_fraction,
				     $rng_seed );
	my $alignment_overlap =
	  $overlap_obj->overlap_fasta_string('');
	my $overlap_length = $overlap_obj->get_overlap_length();
	#	print "id: $qid  overlap_length: $overlap_length  n_sequences: ", $overlap_obj->get_n_sequences(), "\n";
	$string_to_print .= "$overlap_length ";
	if ( $overlap_length >= $min_overlap_length and $overlap_length <= $max_overlap_length and $overlap_obj->get_n_sequences() < $max_n_sequences) { # sufficient overlap
	  my %description_lnL = ();

	  ################## run FastTree  ######################################## 
	  # standard way:
	  my $fasttree_command = "FastTree -wag -gamma -bionj -nosupport $additional_ft_options ";
	  my ($ft_newick, $ft_lnL, $ft_cputime) = run_fasttree($alignment_overlap, $fasttree_command);
	  my $type = "FT";
	  my $description = "$type  $ft_cputime  $ft_newick";
	  $description_lnL{$description} = $ft_lnL;
	  print STDERR "$type $ft_lnL;  ";

	  # using NJ tree as initial tree (rather than bionj done within FT)
	  my ($njft_newick, $njft_lnL, $njft_cputime) = nj_to_ft($alignment_overlap, $alignment_overlap);
	  $type = 'NJ->FT';
	  $description = "$type  $njft_cputime  $njft_newick";
	  $description_lnL{$description} = $njft_lnL;
	  print STDERR "$type $ft_lnL;  ";

	  # using NJ bootstrap trees as initial trees.
	  for my $i_bs (1..$n_bs) {
	    my $bs_alignment_overlap =
	      $overlap_obj->bootstrap_overlap_fasta_string('');
	    my $type = 'BS' . $i_bs . 'NJ->FT ';
	    my ($bsnjft_newick, $bsnjft_lnL, $bsnjft_cputime) = nj_to_ft($bs_alignment_overlap, $alignment_overlap);
	    $description_lnL{"$type  $bsnjft_cputime  $bsnjft_newick"} = $bsnjft_lnL;
	    print STDERR "$type $bsnjft_lnL;  ";
	  }			# end of bootstraps loop
	  print STDERR "\n";

	  # sort by FT likelihood 
	  my @skeys = sort {$description_lnL{$b} <=> $description_lnL{$a} } keys %description_lnL; #Sort by FastTree likelihoods (best first)
	  print "Id $qid  \n";
	  my $best_ft_lnL = $description_lnL{$skeys[0]};
	  my $i = 1;
	  for (@skeys) {
	    my ($descript, $ft_cputime, $newick) = split(" ",$_);
	    my $ft_lnL = $description_lnL{$_};
	    # phyml, if requested.
	    my ($phymlobj, $phymlnewick, $phyml_lnL, $phyml_cput) = ($do_phyml)? 
	      run_phyml($alignment_overlap, $newick, $phyml_opt) : 
		(undef, '()', 0, 0);
	    my $delta_lnL = $best_ft_lnL - $ft_lnL;
	    printf("%3i %28s %10s  %11.2f  %5.3f  %8.6f  %5.2f  %12.2f  %5.1f \n",$i, $qid, $descript, $ft_lnL, $delta_lnL, $delta_lnL/abs($best_ft_lnL), $ft_cputime, $phyml_lnL, $phyml_cput);
	    if ($i == 1) {
	      print $fh_out ("Id $qid  $descript  $ft_lnL  $ft_cputime  $phyml_lnL $phyml_cput \n");
	      print $fh_out "$descript  $newick \n";
	    }
	    $i++;
	  }			# loop over sorted trees
	}
      }
      $state = 'idline';
    } else {			# not an empty line - append to $fasta
      $fasta .= $_;
    }
  }
}


sub nj_to_ft{  # get a tree using NJ and use as init tree for FastTree
  # usually want bootstrap alignment as input to NJ, then real alignment as input to FastTree
  my $nj_input_alignment = shift; #
  my $ft_input_alignment = shift;
  my $nj_newick;
  my ($bs_lnL, $bs_cputime);
  if (1) {			# use NJ for the bs trees ...
    $nj_newick = run_quicktree($nj_input_alignment);
  } else { # can use FastTree (just min. evolution, no ML) for bs trees ...
    ($nj_newick, $bs_lnL, $bs_cputime) = run_fasttree($nj_input_alignment,"FastTree -wag -gamma -bionj -nosupport -noml" );
  }
  my $fasttree_command = "FastTree -wag -gamma -bionj -nosupport ";
  my ($ft_newick, $ft_lnL, $ft_cpu_time) = 
    run_fasttree($ft_input_alignment, $fasttree_command, $nj_newick);
  return ($ft_newick, $ft_lnL, $ft_cpu_time);
}
