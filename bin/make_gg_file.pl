#!/usr/bin/perl -w
use strict;

# the idea here is to make a 'genome-gene association file'
# of the sort used by orthomcl and some of these scripts (ortholog_support.pl ...)
# It takes a set of fasta files, each of which has sequences from
# 1 genome. It finds the genome from the file name, and then 
# writes a file in the gg format.
# At present, assumes a file name of form  Agenome-pep.fasta
# where Agenome is a string which identifies a genome.
# e.g. Alyrata-pep.fasta, Stuberosum-pep.fasta
# Usage example:   make_gg_file.pl > 21species.gg  (no arguments, just reads every file ending in .fasta)
# If your filenames are not like this, you will need to either
# rename them or change the regular expressions used to extract
# the part of the  file name which indicates the genome.
# You may have to change (add to) the following hash which maps
# these bits of file names to species names of form "Genus_species"
# The Genus_species format is not really required; you can use anything as long
# as it is unique. I just wanted something more obviously unambiguous than e.g. 'Brapa' or 'Gmax'
my %spname1_spname2 = (
		       'Alyrata' => 'Arabidopsis_lyrata',
		       'Athaliana' => 'Arabidopsis_thaliana',
		       'Bdistachyon' => 'Brachypodium_distachyon', # grass
		       'Brapa' => 'Brassica_rapa', # turnip
		       'Cpapaya' => 'Carica_papaya', # papaya
		       'Creinhardtii' => 'Chlamydomonas_reinhardtii', # unicellular green alga with flagella
		       'Crubella' => 'Capsella_rubella', # arabidopsis relative
		       'Csativus' => 'Cucumis_sativus', # cucumber
		       'Csinensis' => 'Citrus_sinensis', # orange
		       'Gmax' => 'Glycine_max', # soy
		       'Ljaponicus' => 'Lotus_japonicus',
		       'Mtruncatula' => 'Medicago_truncatula',
		       'Olucimarinus' => 'Ostreococcus_lucimarinus',
		       'Otauri' => 'Ostreococcus_tauri', # unicellular green alga
		       'Osativa' => 'Oryza_sativa', # rice
		       'Ppatens' => 'Physcomitrella_patens', # a moss
		       'Ptrichocarpa' => 'Populus_trichocarpa', # black cottonwood
 		       'Rcommunis' => 'Ricinus_communis', # castorbean
		       'Sbicolor' => 'Sorghum_bicolor', # sorghum
		       'Slycopersicum' => 'Solanum_lycopersicum', # tomato
		       'Smoellendorffii' => 'Selaginella_moellendorffii', # spikemoss
		       'Stuberosum' => 'Solanum_tuberosum', # potato
		       'Thalophila' => 'Thellungiella_halophila', # salt-tolerant relative of arabidopsis
		       'Vcarteri' => 'Volvox_carteri', # green alga; forms spherical colonies.
		       'Vvinifera' => 'Vitis_vinifera', # grape
		       'Zmays' => 'Zea_mays', # maize
		       'Beet' => 'Beta_vulgaris', # sugar beet.
		       'Amborella' => 'Amborella_trichopoda' # basal angiosperm.
);

my %species_sequences_hash = ();
my $fastafile = shift;
my @files;
if(defined $fastafile){
  @files = ($fastafile);
}else{
@files = `ls *.fasta`;
}
map(chomp, @files);
# print "[", join("]\n[", @files), "]\n";
#exit;

foreach my $filename (@files){
	my $species_name = '';
	if($filename =~ /^(.*)[-]pep[.]fasta/){
		$species_name = $1;
		if(exists $spname1_spname2{$species_name}){
		  $species_name =  $spname1_spname2{$species_name}; 
		}
	}else{
		warn "Expected filename of form speciesname-pep.fasta; got: $filename .\n";
		$species_name = '';
		next;
	}
	my @idlines = split("\n", `grep -P '^>' $filename`);
	foreach my $idline (@idlines){
		$idline =~ /^>(\S+)/;
		my $sequence_id = $1;
		$species_sequences_hash{$species_name} .= "$sequence_id ";
#	print scalar keys %species_sequences_hash, "  $species_name   $sequence_id \n";
	}
}
foreach (keys %species_sequences_hash){
	print "$_: ", $species_sequences_hash{$_}, "\n";
}

