#!/usr/bin/perl -w
use strict;

# put taxon names into newick file after each id
# like, e.g.:   AT1G26530.1[species=arabidopsis]
# use lib '/home/tomfy/Orthologger/lib/';
use CXGN::Phylo::IdTaxonMap;

my $id_taxon_map = CXGN::Phylo::IdTaxonMap->new();
my $pattern = shift || '*.newick';

my @inputfiles = split( " ", `ls $pattern` );

foreach my $input_file (@inputfiles) {
    my $output_file = $input_file;
    $output_file =~ s/[.]newick/_taxa.newick/;
    open my $fhin, "<$input_file";
    my $newick = <$fhin>;
    my %map    = %{ $id_taxon_map->get_map() };

    foreach my $id_regex ( keys %map ) {
        my $taxon = $map{$id_regex};
        $id_regex =~ s/^\^//;    # remove initial '^'
        my $tax_str = "[species=$taxon]";
        if ( $newick =~ / ( [\(,] ) \s* ( $id_regex [^,\):]* ) ( [\:] ) /x ) { # ,\):] ) /x ) { 
            print "id_regex: $id_regex ;  1: [$1]  2: [$2]  3: [$3]\n";
        }
        $newick =~ s/ ( [(,] ) \s* ( $id_regex [^,\):]* ) \s* ( [,):] ) /$1$2$tax_str$3/xg;
    }
    open my $fhout, ">$output_file";
    print $fhout "$newick";
    close $fhout;
}
