#!/usr/bin/perl -w
use strict;


my $predefined_taxon_groups =
  { # hashref. keys are names of predef taxon groups; values are hashrefs (keys taxa, values 1)
   '4nonangiosperms' => {
			 'Chlamydomonas_reinhardtii'  => 1,
			 'Physcomitrella_patens'      => 1, 
			 'Selaginella_moellendorffii' => 1, 
		#	 'Pinus_taeda'                => 1, # loblolly pine
			 'Picea_abies' => 1, # norway spruce
			},
   '8basals' => { # which branch off before the monocot-dicot split; Amborella & non-angiosperms.
		'Ostreococcus_tauri' => 1,
		'Ostreococcus_lucimarinus' => 1,
		'Volvox_carteri' => 1,
		'Chlamydomonas_reinhardtii'  => 1,
		'Physcomitrella_patens'      => 1, 
		'Selaginella_moellendorffii' => 1, 
		'Picea_abies' => 1, # norway spruce
		'Amborella_trichopoda' => 1
},
   '8monocots' => {
		   'Phoenix_dactylifera'     => 1, # date palm
		   'Setaria_italica'         => 1, # foxtail millet
		   'Triticum_aestivum'       => 1, # wheat
		   'Hordeam_vulgare'         => 1, # barley
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1 # rice
		  },
   '7monocots' => { # These are the monocots in the 40-species analysis, April-May 2014
		   'Phoenix_dactylifera'     => 1, # date palm
		   'Musa_acuminata' => 1,	   # banana
	      #	   'Setaria_italica'         => 1, # foxtail millet
		   'Triticum_aestivum'       => 1, # wheat
      	   #	   'Hordeum_vulgare'         => 1, # barley
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1 # rice
		  },
   '4monocots' => {
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1 # rice
		  },
  '6monocots' => {
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1, # rice
		  'Musa_acuminata' => 1, # banana
		  'Phoenix_dactylifera' => 1, # date palm
		  },
   '7dicots' => {
		 'Solanum_lycopersicum' => 1, # tomato
		 'Solanum_tuberosum'    => 1, # potato
		 'Vitis_vinifera'       => 1, # grape
		 'Glycine_max'          => 1, # soy
		 'Populus_trichocarpa'  => 1, # poplar
		 'Ricinus_communis'     => 1, # castor
		 'Cucumis_sativus'      => 1  # cucumber
		},
   '8dicots' => {
			     'Solanum_lycopersicum' => 1, # tomato
			     'Solanum_tuberosum'    => 1, # potato
			     'Vitis_vinifera'       => 1, # grape
			     'Glycine_max'          => 1, # soy
			     'Populus_trichocarpa'  => 1, # poplar
			     'Ricinus_communis'     => 1, # castor
			     'Cucumis_sativus'      => 1, # cucumber
			     'Carica_papaya'        => 1  # papaya
			    },
   '5brassicas' => {
		    Brassica_rapa           => 1, # turnip
		    Arabidopsis_thaliana    => 1,
		    Arabidopsis_lyrata      => 1,
		    Thellungiella_halophila => 1,
		    Capsella_rubella        => 1
		   },
   '6negatives' =>  {
		     Brassica_rapa           => 1, # turnip
		     Arabidopsis_thaliana    => 1,
		     Arabidopsis_lyrata      => 1,
		     Thellungiella_halophila => 1,
		     Capsella_rubella        => 1,
		     Beta_vulgaris => 1,
		    },

   '9negatives' =>  {
		     Brassica_rapa           => 1, # turnip
		     Arabidopsis_thaliana    => 1,
		     Arabidopsis_lyrata      => 1,
		     Thellungiella_halophila => 1,
		     Capsella_rubella        => 1,
		     Beta_vulgaris => 1,
		     Nelumbo_nucifera => 1,
		     Utricularia_gibba => 1,
		 'Tarenaya_hassleriana' => 1,
		    },
'12pdicots' => {
		 'Solanum_lycopersicum' => 1, # tomato
		 #	'Solanum_tuberosum'    => 1, # potato
		 'Vitis_vinifera'       => 1,	# grape
		 'Glycine_max'          => 1,	# soy
		 'Populus_trichocarpa'  => 1,	# poplar
		 'Ricinus_communis'     => 1,	# castor
		 'Cucumis_sativus'      => 1,	# cucumber
		 'Aquilegia_coerulea' => 1,	# columbine
		 'Mimulus_guttatus' => 1,	# monkeyflower
		 'Theobroma_cacao' => 1,
		 'Carica_papaya' => 1,

		 'Lupinus_angustifolius' => 1,
		 'Lotus_japonicus' => 1,
		 #	       'Eucalyptus_grandis' => 1,
		 #	       'Manihot_esculenta' => 1,
		},


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
 '9_monocots' => { # These are the monocots in the 50-species analysis Sept. 2014
		       'Panicum_virgatum' => 1, # switchgrass
		       'Phyllostachys_heterocycla' => 1, # bamboo, AM ??
		   'Phoenix_dactylifera'     => 1, # date palm
		   'Musa_acuminata' => 1,	   # banana
		   'Zea_mays'                => 1, # maize
		   'Brachypodium_distachyon' => 1,
		   'Sorghum_bicolor'         => 1,
		   'Oryza_sativa'            => 1, # rice
		    'Spirodela_polyrhiza' => 1, # duckweed - AM negative monocot
    #	   'Setaria_italica'         => 1, # foxtail millet
	#	   'Triticum_aestivum'       => 1, # wheat
      	   #	   'Hordeum_vulgare'         => 1, # barley
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
'7_basals' => { # which branch off before the monocot-dicot split; Amborella & non-angiosperms.
		'Ostreococcus_tauri' => 1,
		'Ostreococcus_lucimarinus' => 1,
		'Volvox_carteri' => 1,
		'Chlamydomonas_reinhardtii'  => 1,
		'Physcomitrella_patens'      => 1, 
		'Selaginella_moellendorffii' => 1, 
		'Picea_abies' => 1, # norway spruce
	#	'Amborella_trichopoda' => 1
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


   # 24 sp for C4 analysis: 
   '9_C3_dicots' => {
		     'Aquilegia_coerulea' => 1,	  # columbine
		     'Solanum_lycopersicum' => 1, # tomato
		     'Vitis_vinifera'       => 1, # grape
		     'Medicago_truncatula' => 1,
		     'Ricinus_communis'     => 1, # castor
		     'Cucumis_sativus'      => 1, # cucumber
		     Arabidopsis_thaliana    => 1,
		     Beta_vulgaris => 1,
		     'Tarenaya_hassleriana' => 1,
		    }, 
   '6_C3_monocots' => {
		       'Brachypodium_distachyon' => 1,
		       'Oryza_sativa' => 1,
		       'Phyllostachys_heterocycla' => 1,
		       'Musa_acuminata' => 1,
		       'Phoenix_dactylifera' => 1,
		       'Spirodela_polyrhiza' => 1,
		      },
 '7_C3_monocots' => {
		       'Brachypodium_distachyon' => 1,
		       'Oryza_sativa' => 1,
		       'Phyllostachys_heterocycla' => 1,
		       'Musa_acuminata' => 1,
		       'Phoenix_dactylifera' => 1,
		       'Spirodela_polyrhiza' => 1,
		     'Hordeum_vulgare' => 1,
		      },
   '4_C4_monocots' => {
		       'Sorghum_bicolor' => 1,
		       'Zea_mays' => 1,
		       'Panicum_virgatum' => 1,
		       'Setaria_italica' => 1,
		      },
   '4_basals' => {
		  'Amborella_trichopoda' => 1,
		  'Picea_abies' => 1,
		  'Selaginella_moellendorffii' => 1,
		  'Physcomitrella_patens' => 1,
		 },
   '19_non_C4s' => {
		    'Amborella_trichopoda' => 1,  # 4 basals
		    'Picea_abies' => 1,
		    'Selaginella_moellendorffii' => 1,
		    'Physcomitrella_patens' => 1,

		    'Aquilegia_coerulea' => 1,	 # columbine   # 9 C3 dicots
		    'Solanum_lycopersicum' => 1, # tomato
		    'Vitis_vinifera'       => 1, # grape
		    'Medicago_truncatula' => 1,
		    'Ricinus_communis'     => 1, # castor
		    'Cucumis_sativus'      => 1, # cucumber
		    Arabidopsis_thaliana    => 1,
		    Beta_vulgaris => 1,
		    'Tarenaya_hassleriana' => 1,

		    'Brachypodium_distachyon' => 1, # 6 C3 monocots
		    'Oryza_sativa' => 1,
		    'Phyllostachys_heterocycla' => 1,
		    'Musa_acuminata' => 1,
		    'Phoenix_dactylifera' => 1,
		    'Spirodela_polyrhiza' => 1,
		   }

  };

my %species_group = ();
my @speciesgroups = ( '8_basals', '9_monocots', '23_AMp_dicots', '11_AMnegatives' );

for my $spgr (@speciesgroups){
  my $species = $predefined_taxon_groups->{$spgr};
  for my $sp (keys %$species){
 #   print "$sp $spgr \n";
    $species_group{$sp} = $spgr;
 #   print STDERR "[$sp]  [$spgr] \n";
  }
}
#exit;

while (<>) {
  if (/^Id \s+ (\S+) \s+ family[.] \s+ fam_size: \s+ (\S+) \s+ (\S+) /x) {
my %speciesgroup_count = ();
for(@speciesgroups){ $speciesgroup_count{$_} = 0; }
# '8basals' => 0, '6monocots' => 0, '12pdicots' => 0, '9negatives' => 0);
    my ($qid, $famsize, $species_string) = ($1, $2, $3);
    print "$qid   $famsize   "; #  $species_string \n";
    my @species_present = split(",", $species_string);
    for(@species_present){
      next if(/Medicago/);
   #   print "$_ ", $species_group{$_}, "\n";
      my $group = $species_group{$_};
      $speciesgroup_count{$group}++;
    }
    for(@speciesgroups){
      print "   ", $speciesgroup_count{$_}, "   ";
    }
    print "\n";
  }
}
