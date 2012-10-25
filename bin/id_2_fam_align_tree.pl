#!/usr/bin/perl -w
use strict;
use warnings;
no lib '/home/tomfy/bin';
no lib '/home/tomfy/Orthologger/bin';
no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ($bindir, $libdir);
BEGIN{ $bindir = dirname(abs_path(__FILE__)); # this has to go in Begin block so happens at compile time
$libdir = $bindir . '/../lib';
}
use lib $libdir;

use Getopt::Long;

use CXGN::Phylo::Parser;
use CXGN::Phylo::BasicTree;
use CXGN::Phylo::File;
use CXGN::Phylo::Species_name_map;

my $default_fasta_file_path   = "/home/tomfy/MHarrison/fall2011families/all.fa";
my $default_cluster_file_path = "/home/tomfy/MHarrison/fall2011families/all_orthomcl.out";

# Defaults
my $id_list_or_file = undef;
my $all_fasta_filename = $default_fasta_file_path; # path to file with sequences in fasta format
my $cluster_filename = $default_cluster_file_path; # path to cluster file (1 line has ids in 1 cluster)
my $do_ortholog_support = 0; # shift || 0; # otherwise just get tree with FastTree
#my $bootstraps_param; #    = shift || '10,10'; # bootstrap params to pass to ortholog_support.pl
my $do_taxonify         = 1;
my $reroot_method = 'none';
my $prune_threshold = 3;

# Process long cl options
GetOptions('ids=s' => \$id_list_or_file,
	   'sequencefile=s' => \$all_fasta_filename,
	   'clusterfile=s' => \$cluster_filename,
	  'reroot=s' => \$reroot_method,
	  'prune_threshold=i' => \$prune_threshold,
	  );
print STDERR "ids: $id_list_or_file\n",
  "reroot method: $reroot_method\n",
  "prune_threshold: $prune_threshold.\n";
           # 'seed=i' => \$seed,
           # 'nongap_fraction=f' => \$nongap_fraction,
           # 'chunk_size=i' => \$chunk_size,
           # 'n_temperatures=i' => \$n_temperatures,
           # 'delta_temperature=f' => \$delta_temperature,
           # 'sample_freq=i' => \$sample_freq,
           # 'n_runs=i' => \$n_runs,
           # 'burn-in_fraction=f' => \$burnin_fraction,
           # 'converged_chunks_required=i' => \$converged_chunks_required,
           # 'ESS_min=i' => \$modelparam_min_ok_ESS);
 ###################
# print "id_list_or_file: [$id_list_or_file]\n";

die "id_2_fam_align_tree.pl requires 3 arguments:\n"
  . " id_or_file (either a seq id or a file containing seq ids 1/line in first column).\n"
  . " x.fasta (sequences, must include id or ids specified by id_or_file).\n"
  . " clusterfile (file w 1 cluster per line).\n"
  if ( !defined $id_list_or_file or !-f $all_fasta_filename or !-f $cluster_filename );
# print "$all_fasta_filename  $cluster_filename \n";

my $id_file = '';
if ( -f $id_list_or_file ) {
  $id_file = $id_list_or_file;
  print STDERR "Getting ids from file: $id_file.\n";
} else {		  # argument is not file, should just be an id
  my $ids = $id_list_or_file;
  print STDERR "Getting ids from list:[$ids]. Sequence, cluster files: $all_fasta_filename  $cluster_filename \n";
$ids =~ s/^\s+//;
$ids =~ s/\s+$//;
$ids =~ s/\s+/\n/g;
$ids = "'" . $ids . "'";
  system  "echo $ids > tmp_idfile";
  $id_file = 'tmp_idfile';
}

my $clusters_ids_2_fasta_cl = $bindir . '/' . 
  "clusters_ids_2_fasta.pl $all_fasta_filename $cluster_filename $id_file";
my $fasta0_filename         = `$clusters_ids_2_fasta_cl`;

my $seqid0;
my @fasta_files = split( "\n", $fasta0_filename );

foreach my $fasta_file (@fasta_files) {
  my @cols = split(" ", $fasta_file);
  $fasta_file = shift @cols;    # just use first col
  $seqid0 = shift @cols;	# could be more than this one.
  my $muscle_input_filename = $fasta_file;
  if ($do_taxonify) {
    my $taxonify_cl = $bindir . '/' . 
      "taxonify_fasta.pl $fasta_file";
    $muscle_input_filename = `$taxonify_cl`;
  }

  #   $muscle_input_filename =~ s/\s*$//;                # remove terminal whitespace

  #    print "rmuscle input file name: [$muscle_input_filename] \n";
  if ( $muscle_input_filename =~ /^\s*(\S+).*$/) {
    $muscle_input_filename = $1;
  } else {
    die "Filename for input to muscle is all whitespace: [$muscle_input_filename].\n";
  }
  my $align_filename = $muscle_input_filename;
  $align_filename =~ s/[.]fasta/_align.fasta/;
  my $muscle_cl  = "muscle -in $muscle_input_filename -out $align_filename";
  my $muscle_out = `muscle -in $muscle_input_filename -out $align_filename`;

  print "Done with alignment using muscle. alignment filename: [$align_filename]\n";

  if ($do_ortholog_support) {

    # print `~/Orthologger/bin/ortholog_support.pl -i $align_filename -T ML -N $bootstraps_param`;
  } else {
    my $FT_cl = "FastTree -wag -gamma -bionj $align_filename 2> $align_filename.out";
    my $FT_out_newick = `$FT_cl`;        

    my $parser = CXGN::Phylo::Parse_newick->new($FT_out_newick, 1);
    my $tree = $parser->parse(CXGN::Phylo::BasicTree->new());

    $tree->get_root()->recursive_implicit_names();
    $tree->get_root()->recursive_implicit_species();
    # if species is empty, (because no [species=x] in newick), get it from the name...
    $tree->set_missing_species_from_names();

    $tree->set_show_standard_species(1);
    $tree->set_species_standardizer(CXGN::Phylo::Species_name_map->new() );

    # rerooting
#    my $opt_r = 'midpoint';
    my $custom_species_tree_newick = undef; # needed only for mindl rerooting. and reroot has default tree for this.
    $tree = reroot($tree, $reroot_method, $custom_species_tree_newick);
    my $whole_tree_newick = $tree->generate_newick();
    print "whole rerooted tree: $whole_tree_newick \n";
    open my $fh_whole_newick, ">$align_filename.newick";
    print $fh_whole_newick  "$whole_tree_newick\n";
    close $fh_whole_newick;

    # pruning
    my $outer_species = {'Selaginella_moellendorffii' => 1,
			 'Zea_mays' => 1, 'Brachypodium_distachyon' => 1, 
			 'Sorghum_bicolor' => 1, 'Oryza_sativa' => 1
			};
    my $pruned_tree_root = $tree->min_clade($seqid0, 3, $outer_species);

    my $pruned_tree_newick = $pruned_tree_root->recursive_generate_newick();
    print "pruned tree: $pruned_tree_newick \n";
    my $pruned_newick_filename = $align_filename . "_pruned.newick";
    open my $fh_pruned_newick, ">$pruned_newick_filename";
    print $fh_pruned_newick  "$pruned_tree_newick\n";
    close $fh_pruned_newick;
  }
}


sub reroot{
  my $tree = shift;
  my $reroot_method = shift || 'midpoint';
  my $species_tree_arg = shift; # used for mindl rerooting. 
  #print "reroot method $reroot_method,  sp tree file: $species_tree_arg \n";

  my ($new_root, $dist_above) = (undef, undef);
  if($reroot_method eq 'none'){
    return $tree;
  } elsif ($reroot_method eq 'midpoint') {
    ($new_root, $dist_above) = $tree->find_midpoint();
  } elsif ($reroot_method eq 'minvar') {
    ($new_root, $dist_above) = $tree->min_leaf_dist_variance_point();
  } elsif ($reroot_method eq 'mindl') {

    my $species_tree_newick;
    if (defined $species_tree_arg) {
      my $species_tree_file = CXGN::Phylo::File->new($species_tree_arg);
      $species_tree_newick = $species_tree_file->get_tree_string();
    } else {
      $species_tree_newick = 
	"(Selaginella[species=Selaginella]:1,(((sorghum[species=Sorghum_bicolor]:1,maize[species=Zea_mays]:1):1,(rice[species=Oryza_sativa]:1,brachypodium[species=Brachypodium_distachyon]:1):1):1,(tomato[species=Solanum_lycopersicum]:1,(grape[species=Vitis_vinifera]:1,((papaya[species=Carica_papaya]:1,arabidopsis[species=Arabidopsis_thaliana]:1):1,((soy[species=Glycine_max]:1,medicago[species=Medicago_truncatula]:1):1,(castorbean[species=Ricinus_communis]:1,Poplar[species=Populus_trichocarpa]:1):1):1):1):1):1):1)";
    }
    my $species_tree_parser = CXGN::Phylo::Parse_newick->new($species_tree_newick, 1);
    my $species_tree = $species_tree_parser->parse(CXGN::Phylo::BasicTree->new());

    $species_tree->set_missing_species_from_names(); # get species from name if species undef
    $species_tree->impose_branch_length_minimum();
    $species_tree->collapse_tree();
    $species_tree->get_root()->recursive_implicit_names();
    $species_tree->get_root()->recursive_implicit_species();

    my $spec_bit_hash = $tree->get_species_bithash($species_tree);

    $tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);
    $species_tree->get_root()->recursive_set_implicit_species_bits($spec_bit_hash);


    ($new_root, $dist_above) = $tree->find_mindl_node($species_tree);
    #  print "after find_mindl. new root node name: ", $new_root->get_name(), " $dist_above\n";
  } else {
    warn "reroot option $reroot_method unknown. No rerooting performed.\n";
  }
  if (defined $new_root) {
    $tree->reset_root_to_point_on_branch($new_root, $dist_above);
  } else {
    warn "In reroot. \$new_root undefined, can't reroot. Tree unchanged.\n";
  }
  return $tree;
}


