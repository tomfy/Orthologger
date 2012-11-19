#!/usr/bin/perl

=head1 NAME

orthologger.pl

=head1 SYNOPSYS
orthologger.pl [ -f -d -x ] -i <newick formatted file> -s <species tree newick file>  
-q <query species> -r <reroot method>
=head1 DESCRIPTION

orthologger.pl determines ortholog relationships from newick formatted trees as input.

=head2 Methodology

Orthologs are found in a gene tree by comparing it with a species tree. Nodes of the gene tree are identified as speciations or duplications, 
and two genes are orthologous if their most recent ancestral node in the gene tree is a speciation. This is dependent on the rooting of the tree. Several rooting methods are options. The recommended rooting method is using the urec program of Gorecki and Tiuryn, which models
the observed gene tree as being the product of speciation, gene duplication and gene loss; for each possible choice of root urec find 
the number of duplication and losses required to give the observed gene tree, and the rooting point is chosen so as to minimize the 
(sum of) duplications and losses. (It is possible to weight the duplications and losses differently in the cost minimized, but the numbers
of duplications and losses is strongly positively correlated so this hardly matters.)

=head2 Specifying the species and the species tree

The determine_species() function in this script maps identifiers to species for the tree supplied with the -i option. This may have to be adapted for different purposes. The species tree needs to be specified in a newick expression, with the nodes using the same species identifiers that the determine_species() function returns. Species tree branch lengths can be set to anything; only the topology of the species tree is used.


=head1 OPTIONS

=over 5

=item -i 

specify a newick input file.

=item -s

the filename of a file containing a Newick formatted tree describing the expected relationships between the species that occur in the sequence tree. Branch lengths in this tree are ignored, but should be !=0.

=item -f

fold subtrees composed of one species into a single node.

=item -d

debug mode. Print a lot of gobbledegook and print results while calculating.

=item -r

Reset root. E.g. -r minvar resets root to point in tree such that variance of root-leaf distances is minimized. Other choices are "minmax", "maxmin". Default is no rerooting.

=item -x

Expansion. Find subtrees that have all leaves the same species, and report number of such leaves.

=back

=head1 VERSION
0.6 2011-12-2
    o Now uses module Orthologger.pm 
0.4 2008-10-3 
    o resetting root in yet a better way (so as to minimize duplications, losses required to reconcile with species tree, using "urec" program of gorecki and tiuryn)
    o finding all orthologs according to definition (two leaves are orthologs if node where lineages come together in gene tree is speciation).
0.3 2007-11-20
    o resetting root in better way [Tom]
    o ortholog group based on leaf count and leaf species count

0.2 2006-08-16
    o added -f option [Scott]

0.1 2005-11-12

=head1 AUTHOR(S)

 Lukas Mueller (lam87@cornell.edu)
 Scott Bland (sbland19@gmail.com)
 Tom York (tly2@cornell.edu)

=head1 LICENSE

Copyright (c) 2002-2006 Lukas Mueller, the Sol Genomics Network and Cornell University.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

=cut
use strict;
use Getopt::Std;

use lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use CXGN::Phylo::File;
use CXGN::Phylo::Node;
use CXGN::Phylo::Tree;
use CXGN::Phylo::Parser;

use lib '/home/tomfy/MHarrisonProject/lib';
use Orthologger;
use vars qw($opt_i $opt_s $opt_d $opt_r $opt_x $opt_q $opt_N);

# get options
getopts("i:s:dr:xq:N:");

if (!$opt_i or !$opt_s) { # gene tree newick file and species tree newick file must both be specified on command line.
  print "Usage: ortholog.pl -i inputfile  -s species_tree [-f -d -x -r <reroot_type; mindl|minvar|maxmin|minmax|midpoint (default=none)>] -q <query species> -g \n\n"; exit();
}

my $species_tree_file_name = $opt_s;
my $reroot_method = ($opt_r)? $opt_r: 'none'; # e.g. "minvar"; by default don't reroot

# print STDERR "in orthologger. reroot method: $reroot_method\n";
my $get_expansion = $opt_x;
my $query_species = $opt_q;

#my $species_name_map = CXGN::Phylo::Species_name_map->new();
# SPECIES TREE:
# now get species tree object from file:
my $species_tree_file = CXGN::Phylo::File->new($species_tree_file_name);
my $species_newick = $species_tree_file->get_tree_string();
# get a parser and parse the newick file. 
my $species_tree = CXGN::Phylo::Parse_newick->new($species_newick)->parse();
if (!$species_tree) {
  die"Species tree. Parse_newick->parse() failed to return a tree object. Newick string: ".$species_newick."\n";
}

foreach my $gene_tree_file_name (`ls $opt_i`) {	# loop over all (gene tree) files specified by $opt_i
  $gene_tree_file_name =~ s/\s//g;

  $opt_N ||= 1;
  foreach (1..$opt_N) {
    print "$_ \n";
    # GENE TREE
    # open the file using CXGN::Phylo::File
    my $gene_tree_file = CXGN::Phylo::File->new($gene_tree_file_name);
    my $gene_tree_newick = $gene_tree_file->get_tree_string(); 
    debug("gene tree newick: $gene_tree_newick \n");
    # get a parser and parse the gene tree newick file; get a tree object. 
    my $gene_tree = CXGN::Phylo::Parse_newick->new($gene_tree_newick)->parse();
    $gene_tree->show_newick_attribute('species'); 
    if (!$gene_tree) {
      die "Gene tree. Parser_newick->parse() failed to return a tree object.\n Newick string: ". $gene_tree_newick . "\n";
    }

    # gene tree object $gene_tree should exist now.

    my $ortholog_obj = Orthologger->new({'gene_tree' => $gene_tree, 'species_tree' => $species_tree, 'reroot_method' => $reroot_method});
# exit;

    print $ortholog_obj->ortholog_result_string(), "\n";
    if ($opt_x) {
      print $ortholog_obj->get_expanded_subtrees($ortholog_obj->get_gene_tree()->get_root(),$gene_tree_file_name );
    }
  }
}

sub debug {
  if ($opt_d) {
    print STDERR join(" ",@_);
  }
}
