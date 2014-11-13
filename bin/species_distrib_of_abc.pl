#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# prints out some information on how many of each species and each
# group of species are present in a family
# input is abc file, gg file, and (optionally) 
# a query id. (If no query id given, gives stats
# for whole abc file.

my $gg_filename = undef;
my $abc_file = undef;
#my $id1 = undef;

my %qid_spcount = ();

my %group_species = (
		     '23_AMp_dicots' => {
					 'Aquilegia_coerulea' => 1, # columbine

					 'Solanum_lycopersicum' => 1, # tomato
					 'Solanum_tuberosum'    => 1, # potato
					 'Mimulus_guttatus' => 1, # monkeyflower
					 'Fraxinus_excelsior' => 1, # Ash
					 'Sesamum_indicum' => 1,

					 'Vitis_vinifera'       => 1, # grape

					 'Glycine_max'          => 1, # soy
					 'Phaseolus_vulgaris' => 1,
					 'Lupinus_angustifolius' => 1,
					 'Lotus_japonicus' => 1,
					 'Medicago_truncatula' => 1,

					 'Populus_trichocarpa'  => 1, # poplar
					 'Ricinus_communis'     => 1, # castor
					 'Cucumis_sativus'      => 1, # cucumber
				  	 'Manihot_esculenta' => 1,
					 'Salix_purpurea' => 1,

					 'Theobroma_cacao' => 1,
					 'Carica_papaya' => 1,
					 'Eucalyptus_grandis' => 1,
					 'Gossypium_raimondii' => 1,
					 'Citrus_clementina' => 1,
					 'Citrus_sinensis' => 1,
					},
		     '8_AMp_monocots' => { # These are the monocots in the 36-species analysis, May 2014
					  'Panicum_virgatum' => 1, # switchgrass
					  'Phyllostachys_heterocycla' => 1, # bamboo, AM ??
					  'Phoenix_dactylifera'     => 1, # date palm
					  'Musa_acuminata' => 1, # banana
					  #	   'Setaria_italica'         => 1, # foxtail millet
					  #	   'Triticum_aestivum'       => 1, # wheat
					  #	   'Hordeum_vulgare'         => 1, # barley
					  'Zea_mays'                => 1, # maize
					  'Brachypodium_distachyon' => 1,
					  'Sorghum_bicolor'         => 1,
					  'Oryza_sativa'            => 1, # rice
		     
					 },
		     '8_basals' => { # which branch off before the monocot-dicot split; Amborella & non-angiosperms.
				    'Ostreococcus_tauri' => 1,
				    'Ostreococcus_lucimarinus' => 1,
				    'Volvox_carteri' => 1,
				    'Chlamydomonas_reinhardtii'  => 1,
				    'Physcomitrella_patens'      => 1, 
				    'Selaginella_moellendorffii' => 1, 
				    'Picea_abies' => 1, # norway spruce
				    'Amborella_trichopoda' => 1
				   },
		     '11_AMnegatives' =>  {
					   Brassica_rapa           => 1, # turnip
					   Arabidopsis_thaliana    => 1,
					   Arabidopsis_lyrata      => 1,
					   Thellungiella_halophila => 1,
					   Capsella_rubella        => 1,
					   Beta_vulgaris => 1,
					   Nelumbo_nucifera => 1,
					   Utricularia_gibba => 1,
					   'Tarenaya_hassleriana' => 1,
					   'Dianthus_caryophyllus' => 1,
					   'Spirodela_polyrhiza' => 1, # duckweed - monocot
					  },
		    );

my %species_group = ();
for my $grp (keys %group_species) {

  my @grp_species = keys %{$group_species{$grp}};
  for (@grp_species) {
    $species_group{$_} = $grp;
  }
}
my %group_seqcount = ();
my %group_speciescount = ();

my @group_list = ('23_AMp_dicots', '8_AMp_monocots', '8_basals', '11_AMnegatives');

my @species = (			# 37
	       # Basal species (8)
	       'Ostreococcus_tauri',
	       'Ostreococcus_lucimarinus',
	       'Volvox_carteri',
	       'Chlamydomonas_reinhardtii',
	       'Physcomitrella_patens',
	       'Selaginella_moellendorffii',
	       'Picea_abies',
	       'Amborella_trichopoda',
	      

	       # monocots: (9;  7 AM+, 1 AM-, 1 AM ??)
	       'Spirodela_polyrhiza', # AM negative
	       'Phoenix_dactylifera',
	       'Musa_acuminata',
	       'Oryza_sativa',
	       'Brachypodium_distachyon',
	       'Zea_mays',
	       'Sorghum_bicolor',
	       'Panicum_virgatum',
	       'Phyllostachys_heterocycla', # AM ??
	       '  ',
	       # AM positive dicots: (23)
	       'Aquilegia_coerulea',
	       'Vitis_vinifera',
	       'Medicago_truncatula',
	       'Glycine_max',
	       'Lupinus_angustifolius',
	       'Lotus_japonicus',
	       'Ricinus_communis',
	       'Cucumis_sativus',
	       'Populus_trichocarpa',
	       'Solanum_lycopersicum',
	       'Mimulus_guttatus',
	       'Theobroma_cacao',
	       'Carica_papaya',
	       'Fraxinus_excelsior', # Ash
	       'Sesamum_indicum',
	       'Solanum_tuberosum',
	       'Phaseolus_vulgaris',
	       'Manihot_esculenta',
	       'Salix_purpurea',
	       'Eucalyptus_grandis',
	       'Gossypium_raimondii',
	       'Citrus_clementina',
	       'Citrus_sinensis',

	       # AM negative dicots: (10)
	       'Arabidopsis_thaliana',
	       'Arabidopsis_lyrata',
	       'Capsella_rubella',
	       'Thellungiella_halophila',
	       'Brassica_rapa',
	       'Tarenaya_hassleriana',
	       'Beta_vulgaris',
	       'Nelumbo_nucifera',
	       'Utricularia_gibba',
	       'Dianthus_caryophyllus',
	      ); 


GetOptions(
	   'gg_file=s'           => \$gg_filename,
	   'abc_file=s'          => \$abc_file,
#	   'id1=s' => \$id1,
	  );

my $gg = store_gene_genome_association_info($gg_filename);
my %species_count = ();
open my $fh, "<", "$abc_file" or die " couldnt open $abc_file for reading. \n";


while (<$fh>) {
  next if(/^\s*$/);
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my $qid = $cols[0];
  # next if(defined $id1 and ($id1 ne $cols[0]));

  my $species = $gg->{$cols[1]};
  if (exists $qid_spcount{$qid}) {
    $qid_spcount{$qid}->{$species}++;
  } else {
    $qid_spcount{$qid} = {$species => 1};
  }
  $species_count{$species}++;
  # print $cols[1], "  ", $species, "\n";
  $group_seqcount{$species_group{$species}}++;
}
# for(keys %qid_spcount){
#   print "AA: $_   ", join(", ", keys %{$qid_spcount{$_}}), "\n";
# }


for my $qid (keys %qid_spcount) {
  my $sp_count = $qid_spcount{$qid};
  my $species_count_string = '';
my $i = 0;
  for my $spec (@species) {
    my $sp_count_x = (exists $sp_count->{$spec})? $sp_count->{$spec} : 0;
    $sp_count_x = 9 if($sp_count_x > 9);
    $species_count_string .= $sp_count_x;
    $i++;
    $species_count_string .= ' ' if($i == 9   or  $i == 17  or $i == 40);
  }
  print "$species_count_string  $qid\n";
}

sub store_gene_genome_association_info {
  my $gg_filename = shift;
  my %gene_genome = ();

  open my $fh, "<", "$gg_filename";
  while (<$fh>) {
    my @cols = split( " ", $_ );
    my $genome = shift @cols;
    $genome =~ s/:$//;		# remove the colon.
    for my $the_id (@cols) {
      $gene_genome{$the_id} = $genome;
    }
  }
  return \%gene_genome;
}
