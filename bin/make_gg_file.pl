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
# these file names to species names of form "Genus_species"
# The Genus_species format is not really required; you can use anything as long
# as it is unique. I just wanted something more obviously unambiguous than e.g. 'Brapa' or 'Gmax'
my %filename_spname = (
'Acoerulea_195_v1.1.protein_primaryTranscriptOnly.fasta' => 'Aquilegia_coerulea',
'Alyrata-pep.fasta' => 'Arabidopsis_lyrata',
'AmTr_v1.0_evm_run27_filter02.prot.fasta' => 'Amborella_trichopoda',
'AshUniproteins.fasta' => 'Fraxinus_excelsior',
'Athaliana-pep.fasta' => 'Arabidopsis_thaliana',
'Bdistachyon_283_v2.1.protein_primaryTranscriptOnly.fasta' => 'Brachypodium_distachyon',
'Brapa-pep.fasta' => 'Brassica_rapa',
'cacao11genes_pub3i.aa.fasta' => 'Theobroma_cacao',
'Cclementina_182_v1.0.protein_primaryTranscriptOnly.fasta' => 'Citrus_clementina',
'Cpapaya_113_ASGPBv0.4.protein_primaryTranscriptOnly.fasta' => 'Carica_papaya',
'Creinhardtii-pep.fasta' => 'Chlamydomonas_reinhardtii',
'Crubella-pep.fasta' => 'Capsella_rubella',
'Csativus-pep.fasta' => 'Cucumis_sativus',
'Csinensis_154_v1.1.protein_primaryTranscriptOnly.fasta' => 'Citrus_sinensis',
'DCA_r1.0_pep.fasta' => 'Dianthus_caryophyllus',
'Egrandis_201_v1.1.protein_primaryTranscriptOnly.fasta' => 'Eucalyptus_grandis',
'Esalsugineum_173_v1.0.protein_primaryTranscriptOnly.fasta' => 'Thellungiella_halophila',
'Gmax_189_protein_primaryTranscriptOnly.fasta' => 'Glycine_max',
'Graimondii_221_v2.1.protein_primaryTranscriptOnly.fasta' => 'Gossypium_raimondii',
'Lj2.5_proteins.fasta' => 'Lotus_japonicus',
'Lupin.pep.fasta' => 'Lupinus_angustifolius',
'Mesculenta-pep.fasta' => 'Manihot_esculenta',
'Mguttatus_256_v2.0.protein_primaryTranscriptOnly.fasta' => 'Mimulus_guttatus',
'Mt4.0v1_GenesProteinSeq_20130731_1800.fasta' => 'Medicago_truncatula',
'Musa_acuminata.MA1.21.pep.all.fasta' => 'Musa_acuminata',
'Olucimarinus-pep.fasta' => 'Ostreococcus_lucimarinus',
'Osativa_204_v7.0.protein_primaryTranscriptOnly.fasta' => 'Oryza_sativa',
'Otauri-pep.fasta' => 'Ostreococcus_tauri',
'Pabies1.0-all-pep.fasta' => 'Picea_abies',
'Pdactylifera-pep.fasta' => 'Phoenix_dactylifera',
'P_heterocycla_v1.0.genemodel.protein.fasta' => 'Phyllostachys_heterocycla',
'Ppatens_251_v3.0.protein_primaryTranscriptOnly.fasta' => 'Physcomitrella_patens',
'Ptrichocarpa_210_v3.0.protein_primaryTranscriptOnly.fasta' => 'Populus_trichocarpa',
'Pvirgatum_273_v1.1.protein_primaryTranscriptOnly.fasta' => 'Panicum_virgatum',
'Pvulgaris_218_v1.0.protein_primaryTranscriptOnly.fasta' => 'Phaseolus_vulgaris',
'Rcommunis-pep.fasta' => 'Ricinus_communis',
'RefBeet.genes.1302.pep.fasta' => 'Beta_vulgaris',
'Sacred_Lotus-pep.fasta' => 'Nelumbo_nucifera',
'Sbicolor_255_v2.1.protein_primaryTranscriptOnly.fasta' => 'Sorghum_bicolor',
'Sesamum_indicum_v1.0.gene.pep.fasta' => 'Sesamum_indicum',
'Slycopersicum-pep.fasta' => 'Solanum_lycopersicum',
'Smoellendorffii_91_v1.0.protein_primaryTranscriptOnly.fasta' => 'Selaginella_moellendorffii',
'Spolyrhiza_290_v2.protein_primaryTranscriptOnly.fasta' => 'Spirodela_polyrhiza',
'Spurpurea_289_v1.0.protein_primaryTranscriptOnly.fasta' => 'Salix_purpurea',
'Stuberosum-pep.fasta' => 'Solanum_tuberosum',
'Thassleriana-pep.fasta' => 'Tarenaya_hassleriana',
'Ugibba-prot.fasta' => 'Utricularia_gibba',
'Vcarteri-pep.fasta' => 'Volvox_carteri',
'Vvinifera_145_Genoscope.12X.protein_primaryTranscriptOnly.fasta' => 'Vitis_vinifera',
'Zmays_284_6a.protein_primaryTranscriptOnly.fasta' => 'Zea_mays',
);

my %C4file_sp = (
		 'Acoerulea-pep.fasta' => 'Aquilegia_coerulea',
		 'Athaliana-pep.fasta' => 'Arabidopsis_thaliana',
		 'Atrichopoda-pep.fasta' => 'Amborella_trichopoda',
		 'Bdistachyon-pep.fasta' => 'Brachypodium_distachyon',
		 'Bvulgaris-pep.fasta' => 'Beta_vulgaris',
		 'Csativus-pep.fasta' => 'Cucumis_sativus',
		 'Macuminata-pep.fasta' => 'Musa_acuminata',
		 'Mtruncatula-v4-pep.fasta' => 'Medicago_truncatula',
		 'Osativa-pep.fasta' => 'Oryza_sativa',
		 'P_heterocycla_v1.0.genemodel.protein.fasta' => 'Phyllostachys_heterocycla', # bamboo
		 'Pabies-pep.fasta' => 'Picea_abies',
		 'Pdactylifera-pep.fasta' => 'Phoenix_dactylifera',
		 'Ppatens-pep.fasta' => 'Physcomitrella_patens',
		 'Pvirgatum_273_v1.1.protein_primaryTranscriptOnly.fasta' => 'Panicum_virgatum',
		 'Rcommunis-pep.fasta' => 'Ricinus_communis',
		 'Sbicolor-pep.fasta' => 'Sorghum_bicolor',
		 'Setaria_italica.JGIv2.0.23.pep.all.fasta' => 'Setaria_italica',
		 'Smoellendorffii-pep.fasta' => 'Selaginella_moellendorffii',
		 'Solyc_ITAG2.4_proteins.fasta' => 'Solanum_lycopersicum',
		 'Spolyrhiza_290_v2.protein_primaryTranscriptOnly.fasta' => 'Spirodela_polyrhiza', 
		 'Thassleriana-pep.fasta' => 'Tarenaya_hassleriana',
		 'Vvinifera-pep.fasta' => 'Vitis_vinifera',
		 'Zmays-pep.fasta' => 'Zea_mays',
		);

while (my ($f,$sp) = each %C4file_sp) {
  if (exists $filename_spname{$f}) {
    print STDERR "YY filename, species: $f $sp \n";
  } else {
    $filename_spname{$f} = $sp;
    print STDERR "XX filename, species: $f $sp \n";
  }
}
;

# 		       # AM negative
# 		       'Alyrata' => 'Arabidopsis_lyrata',
# 		       'Athaliana' => 'Arabidopsis_thaliana',
# 		       'Brapa' => 'Brassica_rapa', # turnip
# 		       'Crubella' => 'Capsella_rubella', # arabidopsis relative
# 		       'Thalophila' => 'Thellungiella_halophila', # salt-tolerant relative of arabidopsis
# 		       'Bvulgaris' => 'Beta_vulgaris', # sugar beet.
# 		       'Nnucifera' => 'Nelumbo_nucifera',
# 		       'Ugibba' => 'Utricularia_gibba', 
# 		       'Thassleriana' => 'Tarenaya_hassleriana',
# 		       'Spolyrhiza' => 'Spirodela_polyrhiza', # monocot

# 		       # AM positive monocots
# 		       'Bdistachyon' => 'Brachypodium_distachyon', # grass
# 		       'Osativa' => 'Oryza_sativa', # rice
# 		       'Sbicolor' => 'Sorghum_bicolor', # sorghum
# 		       'Zmays' => 'Zea_mays',		# maize
# 		       'Macuminata' => 'Musa_acuminata',
# 		       'Pdactylifera' => 'Phoenix_dactylifera',

# 		       # AM positive dicots
# 		       'Acoerulea' => 'Aquilegia_coerulea',
# 		       'Cpapaya' => 'Carica_papaya', # papaya
#         	       'Csativus' => 'Cucumis_sativus',	 # cucumber
# 		       'Gmax' => 'Glycine_max',		 # soy
# 		       'Ljaponicus' => 'Lotus_japonicus',
# 		       #    'Mtruncatula' => 'Medicago_truncatula',
# 		       'MtruncatulaJCVI4.0' => 'Medicago_truncatula',
# 		       'Ptrichocarpa' => 'Populus_trichocarpa', # black cottonwood
#  		       'Rcommunis' => 'Ricinus_communis', # castorbean
# 		       #    'Slycopersicum' => 'Solanum_lycopersicum', # tomato
# 		       'Slycopersicum_ITAG2.4' => 'Solanum_lycopersicum', 
# 		       'Vvinifera' => 'Vitis_vinifera',	    # grape
# 		       'Langustifolius' => 'Lupinus_angustifolius', # lupine
# 		       'Tcacao' => 'Theobroma_cacao', # cacao
# 		       'Mguttatus' => 'Mimulus_guttatus', # monkeyflower

# 		       # basal species (amborella and non-angiosperms
# 		       'Amborella' => 'Amborella_trichopoda' # basal angiosperm.
# 		       'Creinhardtii' => 'Chlamydomonas_reinhardtii', # unicellular green alga with flagella
# 		       'Olucimarinus' => 'Ostreococcus_lucimarinus',
# 		       'Otauri' => 'Ostreococcus_tauri', # unicellular green alga
# 		       'Ppatens' => 'Physcomitrella_patens', # a moss
# 		       'Pabies' => 'Picea_abies',
# 		       'Smoellendorffii' => 'Selaginella_moellendorffii', # spikemoss
# 		       'Vcarteri' => 'Volvox_carteri', # green alga; forms spherical colonies.

# # not in 37 species set
#  'Csinensis' => 'Citrus_sinensis', # orange
#  'Stuberosum' => 'Solanum_tuberosum', # potato
# 		      );

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

foreach my $filename (keys %filename_spname){  # @files) {
  my $species_name = '';
  #	if($filename =~ /^(.*)[-]pep[.]fasta/){
  if (-f $filename) {
    #	$species_name = $1;
    if (exists $filename_spname{$filename}) {
      print STDERR "$filename    $species_name \n";
      $species_name =  $filename_spname{$filename};
    } else {
      warn "Fasta file: $filename has no associated species name.\n";
    }
  } else {
    warn "No file with name $filename .\n";
    $species_name = '';
    next;
  }
  my @idlines = split("\n", `grep -P '^>' $filename`);
  foreach my $idline (@idlines) {
    $idline =~ /^>(\S+)/;
    my $sequence_id = $1;
    $species_sequences_hash{$species_name} .= "$sequence_id ";
    #	print scalar keys %species_sequences_hash, "  $species_name   $sequence_id \n";
  }
}
foreach (keys %species_sequences_hash) {
  print "$_: ", $species_sequences_hash{$_}, "\n";
}

