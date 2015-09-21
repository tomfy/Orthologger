#!/usr/bin/perl -w
use strict;

my $taxon_groups = {
                    '22_AM_dicots' => {
                                       'Aquilegia_coerulea' => 1, # columbine

                                       'Solanum_lycopersicum' => 1, # tomato
                                       'Solanum_tuberosum'    => 1, # potato
                                       'Mimulus_guttatus' => 1, # monkeyflower
                                       'Fraxinus_excelsior' => 1, # Ash
                                       'Sesamum_indicum' => 1,

                                       'Vitis_vinifera'       => 1, # grape

                                       'Glycine_max'          => 1, # soy
                                       'Phaseolus_vulgaris' => 1,
                                       #  'Lupinus_angustifolius' => 1,
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

                    '12_AM_monocots' => { # These are the monocots in the 50-species analysis Sept. 2014
                                         # 8 AM + nonAM Spirodela
                                         'Musa_acuminata' => 1, # banana
                                         'Phoenix_dactylifera'     => 1, # date palm
                                         'Elaeis_guineensis' => 1,

                                         'Phyllostachys_heterocycla' => 1, # bamboo, AM ??
                                         'Zea_mays'                => 1, # maize
                                         'Brachypodium_distachyon' => 1,
                                         'Sorghum_bicolor'         => 1,
                                         'Oryza_sativa'            => 1, # rice
                                         'Panicum_virgatum' => 1, # switchgrass
                                         'Panicum_hallii' => 1,
                                         'Triticum_urartu' => 1,
                                         'Hordeum_vulgare' => 1,
                                        }
                   };

while (<>) {
   if (/^Id\s+(\S+)\s+.*fam_size:\s+(\S+)\s+(\S+)/) {
      my $id = $1;
      my $fam_size = $2;
      my @species = split(",", $3);
      my $AM_dicot_species_count = 0;
      my $AM_monocot_species_count = 0;
      for (@species) {
         $AM_dicot_species_count++ if(exists $taxon_groups->{'22_AM_dicots'}->{$_});
         $AM_monocot_species_count++ if(exists $taxon_groups->{'12_AM_monocots'}->{$_});
      }
 print "$id   $AM_dicot_species_count  $AM_monocot_species_count \n";
   }
  
}

