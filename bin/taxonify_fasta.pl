#!/usr/bin/perl -w
use strict;

# put taxon names into the fasta file id
# e.g. AT1G26530[species=arabidopsis]
# use lib '/home/tomfy/Orthologger/lib';

no lib '/home/tomfy/bin';
no lib '/home/tomfy/Orthologger/bin';
no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use FindBin qw($Bin);
use lib "$Bin/../lib";

use CXGN::Phylo::IdTaxonMap;

my $pattern = shift; # || '*.fasta';
die "No input file specified. Usage example:  taxonify_fasta '*.fasta'" unless (defined $pattern);
my @files = split(" ", `ls $pattern`);

my $id_taxon_map = CXGN::Phylo::IdTaxonMap->new();

foreach my $input_file (@files){

	open my $fhin, "<$input_file";
	my @lines = <$fhin>;
	close $fhin;

	my $out_string = '';
	foreach (@lines) {
		if (/^>/) {
			my $id = $_;
			chomp $id;
			$id = fix_id($id);
			$id =~ s/^>//;
			my $taxon_name = $id_taxon_map->id_to_taxonname($id);
			$id = '>' . $id . "[species=$taxon_name]";
			$out_string .= "$id\n";
		} else {
			$out_string .= $_;
		}
	}

	my $output_file = $input_file;
	$output_file =~ s/[.]fasta/_tax.fasta/;
	if($output_file eq $input_file){
		$output_file .= '_tax';
	}
	open my $fhout, ">$output_file";
	print $fhout $out_string;
	close $fhout;
	print #"Created output file: ", 
		$output_file, "\n";
}



sub fix_id{ # this does a couple of fixes to the id name
#  IMGA|... -> IMGA_... and back to IMGA| after del rest of stuff after |

	my $id_string = shift;
# fixes to $align_string:
	$id_string =~ s/IMGA[|]/IMGA_/g; #pipes in id cause problem; replace '|' with '_'.
		$id_string =~ s/(>[^|]+)[|].*/$1/g; # delete everything from first pipe on.
$id_string =~ s/^>?IMGA_/IMGA|/;
#print STDERR "In taxonify. IMGA=containing id: $id_string \n" if($id_string =~ /IMGA/);

#		my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this. Only needed for 'clearcut' NJ program - so don't need now.
#		$id_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg;

	return $id_string;
}
