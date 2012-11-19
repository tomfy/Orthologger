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

my $using_gg = 0;
my %seqid_species = ();
my $gg_file = shift;
print STDERR join(", ", @files), "  $gg_file \n";
if(defined $gg_file and -f $gg_file){
  $using_gg = 1;
  open my $fh_gg, "<$gg_file";
  while(<$fh_gg>){
    my @cols = split(" ", $_);
    my $species = shift @cols;
    $species =~ s/:$//; # remove final colon if present.
    for (@cols){
      $seqid_species{$_} = $species;
    }
  }
}
print STDERR "Done with storing info from gg file.\n";
my $id_taxon_map = CXGN::Phylo::IdTaxonMap->new();

foreach my $input_file (@files){

	open my $fhin, "<$input_file";
	my @lines = <$fhin>;
	close $fhin;

	my $out_string = '';
	foreach (@lines) {
	  if (/^>/) {
	    /^\s*(\S+)/;
	    my $id = $1;
	    $id = fix_id($id);
	    $id =~ s/^>//;
	    my $taxon_name;
	    if ($using_gg) {
	      if (exists $seqid_species{$id}) {
		$taxon_name = $seqid_species{$id};
	      } else {
		warn "no taxon found for [$id] using gg file.\n";
		$taxon_name = $id_taxon_map->id_to_taxonname($id);
	      }
	    } else {
	      $taxon_name = $id_taxon_map->id_to_taxonname($id);
	    }
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
		$id_string =~ s/^(>?[^|]+)[|].*/$1/g; # delete everything from first pipe on.
$id_string =~ s/^>?IMGA_/IMGA|/;
#print STDERR "In taxonify. IMGA=containing id: $id_string \n" if($id_string =~ /IMGA/);

#		my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this. Only needed for 'clearcut' NJ program - so don't need now.
#		$id_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg;

	return $id_string;
}
