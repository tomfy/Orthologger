#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# prints out some information on how many of each species and each
# group of species are present in a family
# input is abc file (or .m8 file), gg file, and (optionally) 
# a query id. (If no query id given, gives stats
# for whole abc file.

      my $gg_filename = undef;
      my $abc_file = undef;
      my $id1 = undef;

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
         '22_AMp_dicots' => {
            'Aquilegia_coerulea' => 1, # columbine

               'Solanum_lycopersicum' => 1, # tomato
               'Solanum_tuberosum'    => 1, # potato
               'Mimulus_guttatus' => 1, # monkeyflower
               'Fraxinus_excelsior' => 1, # Ash
               'Sesamum_indicum' => 1,

            'Vitis_vinifera'       => 1, # grape

               'Glycine_max'          => 1, # soy
               'Phaseolus_vulgaris' => 1,
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
         '12_AMp_monocots' => { # These are the monocots in the 36-species analysis, May 2014
            'Panicum_virgatum' => 1, # switchgrass
               'Panicum_hallii' => 1,		
            'Phyllostachys_heterocycla' => 1, # bamboo, AM ??

               'Phoenix_dactylifera'     => 1, # date palm
               'Musa_acuminata' => 1, # banana
               'Elaeis_guineensis' => 1,
#	   'Setaria_italica'         => 1, # foxtail millet
#	   'Triticum_aestivum'       => 1, # wheat
            
'Hordeum_vulgare'         => 1, # barley
               'Triticum_urartu' => 1,				  
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
         '13_AMnegatives' =>  {
            Brassica_rapa           => 1, # turnip
               Arabidopsis_thaliana    => 1,
            Arabidopsis_lyrata      => 1,
#  Thellungiella_halophila => 1,
            Eutrema_salsugineum => 1, 		  
            Capsella_rubella        => 1,

            Beta_vulgaris => 1,
            Nelumbo_nucifera => 1,
            Utricularia_gibba => 1,
            'Tarenaya_hassleriana' => 1,
            'Dianthus_caryophyllus' => 1,

            'Spinacia_oleracea' => 1,
            'Spirodela_polyrhiza' => 1, # duckweed - monocot
               'Zostera_marina' => 1, # eelgrass			
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

my @group_list = (
# '23_AMp_dicots', 
      '22_AMp_dicots', 
      '12_AMp_monocots', '8_basals', '13_AMnegatives');

my @species = ();
for my $grp (@group_list){
for my $sp (keys %{$group_species{$grp}}){
   push @species, $sp;
}
}

my @species1 = (			# 37
# Basal species (8)
      'Ostreococcus_tauri',
      'Ostreococcus_lucimarinus',
      'Volvox_carteri',
      'Chlamydomonas_reinhardtii',
      'Physcomitrella_patens',
      'Selaginella_moellendorffii',
      'Picea_abies',
      'Amborella_trichopoda',
      ' ', 

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
      '   ',
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
      'id1=s' => \$id1,
      );

my $gg = store_gene_genome_association_info($gg_filename);
my %species_count = ();
open my $fh, "<", "$abc_file" or die " couldnt open $abc_file for reading. \n";


while (<$fh>) {
   next if(/^\s*$/);
   next if(/^\s*#/);
   my @cols = split(" ", $_);
   next if(defined $id1 and ($id1 ne $cols[0]));
   my $species = $gg->{$cols[1]};
   $species_count{$species}++;
# print $cols[1], "  ", $species, "\n";
   $group_seqcount{$species_group{$species}}++;
}
close $fh;

my $species_content_string = '';
for (@species) {
   my $count = (defined $species_count{$_})? $species_count{$_} : 0;
   my $count_1digit = ($count >= 10)? 9 : $count; 
   $species_content_string .= $count_1digit;
   if(/\S+/){
      printf("%36s   %8i \n", $_, $count);
   }else{
      print  "\n";
   }
   if ($count > 0) {
      my $group = $species_group{$_};
      $group_speciescount{$group}++;
#  print "XXX $_  ", $group_speciescount{$group}, "\n";
   }
}
print "\n";
print "$species_content_string \n";

my $all_groups_seq_count = 0;
my $all_groups_species_count = 0;
for (@group_list) {
   my $seq_count = (defined $group_seqcount{$_})? $group_seqcount{$_}: 0;
   my $species_count = (defined $group_speciescount{$_})? $group_speciescount{$_}: 0;
   printf ("%30s  %8i %8i \n", $_, $seq_count, $species_count);
   $all_groups_seq_count += $seq_count;
   $all_groups_species_count += $species_count;
}

printf ("                        Totals  %8i %8i \n", $all_groups_seq_count, $all_groups_species_count);








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
