#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# read in a tree (newick) format,
# from a newicks file
# and write out a nexus format file for 
# use by figtree, with colors specified
# which depend on species

my $gg_filename = undef;
my $newicks_filename = undef;
my $id = undef;
my $keep_branch_supports = 1; # default: keep branch supports in newick expression
my $grps = 'AM'; # or C4
GetOptions(
	   'gg_file=s'      => \$gg_filename, # not needed if species info present in newicks file.
	   'newicks_file=s' => \$newicks_filename,
	   'id=s' => \$id,
	   'groups=s' => \$grps,
	   
	  );

my $figtree_block = 'begin figtree;' . "\n" . 
  '    set appearance.backgroundColorAttribute="Default";' . "\n" .
  '    set appearance.backgroundColour=#ffffff;' . "\n" .
  '    set appearance.branchColorAttribute="User selection";' . "\n" .
  '    set appearance.branchColorGradient=false;' . "\n" .
  '    set appearance.branchLineWidth=1.0;' . "\n" .
  '    set appearance.branchMinLineWidth=0.0;' . "\n" .
  '    set appearance.branchWidthAttribute="Fixed";' . "\n" .
  '    set appearance.foregroundColour=#000000;' . "\n" .
  '    set appearance.hilightingGradient=false;' . "\n" .
  '    set appearance.selectionColour=#2d3680;' . "\n" .
  '    set branchLabels.colorAttribute="User selection";' . "\n" .
  '    set branchLabels.displayAttribute="bs";' . "\n" .
  '    set branchLabels.fontName="Arial";' . "\n" .
  '    set branchLabels.fontSize=11;' . "\n" .
  '    set branchLabels.fontStyle=0;' . "\n" .
  '    set branchLabels.isShown=true;' . "\n" .
  '    set branchLabels.significantDigits=4;' . "\n" .
  '    set layout.expansion=0;' . "\n" .
  '    set layout.layoutType="RECTANGULAR";' . "\n" .
  '    set layout.zoom=0;' . "\n" .
  '    set legend.attribute="label";' . "\n" .
  '    set legend.fontSize=10.0;' . "\n" .
  '    set legend.isShown=false;' . "\n" .
  '    set legend.significantDigits=4;' . "\n" .
  '    set nodeBars.barWidth=4.0;' . "\n" .
  '    set nodeBars.displayAttribute=null;' . "\n" .
  '    set nodeBars.isShown=false;' . "\n" .
  '    set nodeLabels.colorAttribute="User selection";' . "\n" .
  '    set nodeLabels.displayAttribute="Node ages";' . "\n" .
  '    set nodeLabels.fontName="Arial";' . "\n" .
  '    set nodeLabels.fontSize=11;' . "\n" .
  '    set nodeLabels.fontStyle=0;' . "\n" .
  '    set nodeLabels.isShown=false;' . "\n" .
  '    set nodeLabels.significantDigits=4;' . "\n" .
  '    set nodeShape.colourAttribute="User selection";' . "\n" .
  '    set nodeShape.isShown=false;' . "\n" .
  '    set nodeShape.minSize=10.0;' . "\n" .
  '    set nodeShape.scaleType=Width;' . "\n" .
  '    set nodeShape.shapeType=Circle;' . "\n" .
  '    set nodeShape.size=4.0;' . "\n" .
  '    set nodeShape.sizeAttribute="Fixed";' . "\n" .
  '    set polarLayout.alignTipLabels=false;' . "\n" .
  '    set polarLayout.angularRange=0;' . "\n" .
  '    set polarLayout.rootAngle=0;' . "\n" .
  '    set polarLayout.rootLength=100;' . "\n" .
  '    set polarLayout.showRoot=true;' . "\n" .
  '    set radialLayout.spread=0.0;' . "\n" .
  '    set rectilinearLayout.alignTipLabels=false;' . "\n" .
  '    set rectilinearLayout.curvature=0;' . "\n" .
  '    set rectilinearLayout.rootLength=100;' . "\n" .
  '    set scale.offsetAge=0.0;' . "\n" .
  '    set scale.rootAge=1.0;' . "\n" .
  '    set scale.scaleFactor=1.0;' . "\n" .
  '    set scale.scaleRoot=false;' . "\n" .
  '    set scaleAxis.automaticScale=true;' . "\n" .
  '    set scaleAxis.fontSize=8.0;' . "\n" .
  '    set scaleAxis.isShown=false;' . "\n" .
  '    set scaleAxis.lineWidth=1.0;' . "\n" .
  '    set scaleAxis.majorTicks=1.0;' . "\n" .
  '    set scaleAxis.origin=0.0;' . "\n" .
  '    set scaleAxis.reverseAxis=false;' . "\n" .
  '    set scaleAxis.showGrid=true;' . "\n" .
  '    set scaleBar.automaticScale=true;' . "\n" .
  '    set scaleBar.fontSize=10.0;' . "\n" .
  '    set scaleBar.isShown=true;' . "\n" .
  '    set scaleBar.lineWidth=1.0;' . "\n" .
  '    set scaleBar.scaleRange=0.0;' . "\n" .
  '    set tipLabels.colorAttribute="User selection";' . "\n" .
  '    set tipLabels.displayAttribute="Names";' . "\n" .
  '    set tipLabels.fontName="Arial";' . "\n" .
  '    set tipLabels.fontSize=12;' . "\n" .
  '    set tipLabels.fontStyle=0;' . "\n" .
  '    set tipLabels.isShown=true;' . "\n" .
  '    set tipLabels.significantDigits=4;' . "\n" .
  '    set trees.order=false;' . "\n" .
  '    set trees.orderType="increasing";' . "\n" .
  '    set trees.rooting=true;' . "\n" .
  '    set trees.rootingType="User Selection";' . "\n" .
  '    set trees.transform=false;' . "\n" .
  '    set trees.transformType="cladogram";' . "\n" .
  'end;' . "\n";


my $AM_groups = {
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
	      '8_AMp_monocots' => { # These are the monocots in the 50-species analysis Sept. 2014
				   'Panicum_virgatum' => 1, # switchgrass
				   'Phyllostachys_heterocycla' => 1, # bamboo, AM ??
				   'Phoenix_dactylifera'     => 1, # date palm
				   'Musa_acuminata' => 1, # banana
				   'Zea_mays'                => 1, # maize
				   'Brachypodium_distachyon' => 1,
				   'Sorghum_bicolor'         => 1,
				   'Oryza_sativa'            => 1, # rice
				   #	  
				   #	   '    setaria_italica'         => 1, # foxtail millet
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
	     };

my $C4_groups =  {  # 24 sp for C4 analysis: 
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
   # '6_C3_monocots' => {
   # 		       'Brachypodium_distachyon' => 1,
   # 		       'Oryza_sativa' => 1,
   # 		       'Phyllostachys_heterocycla' => 1,
   # 		       'Musa_acuminata' => 1,
   # 		       'Phoenix_dactylifera' => 1,
   # 		       'Spirodela_polyrhiza' => 1,
   # 		      },
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
   # '19_non_C4s' => {
   # 		    'Amborella_trichopoda' => 1,  # 4 basals
   # 		    'Picea_abies' => 1,
   # 		    'Selaginella_moellendorffii' => 1,
   # 		    'Physcomitrella_patens' => 1,

   # 		    'Aquilegia_coerulea' => 1,	 # columbine   # 9 C3 dicots
   # 		    'Solanum_lycopersicum' => 1, # tomato
   # 		    'Vitis_vinifera'       => 1, # grape
   # 		    'Medicago_truncatula' => 1,
   # 		    'Ricinus_communis'     => 1, # castor
   # 		    'Cucumis_sativus'      => 1, # cucumber
   # 		    Arabidopsis_thaliana    => 1,
   # 		    Beta_vulgaris => 1,
   # 		    'Tarenaya_hassleriana' => 1,

   # 		    'Brachypodium_distachyon' => 1, # 6 C3 monocots
   # 		    'Oryza_sativa' => 1,
   # 		    'Phyllostachys_heterocycla' => 1,
   # 		    'Musa_acuminata' => 1,
   # 		    'Phoenix_dactylifera' => 1,
   # 		    'Spirodela_polyrhiza' => 1,
   # 		   }

  };

my $groups = ($grps eq 'AM')? $AM_groups : $C4_groups;

my $qcolor = '#00ff00';
my %group_color = ('11_AMnegatives' => '#ff0000',
		   '23_AMp_dicots' => '#000000',
		   '8_basals' => '#0000ff',
		   '8_AMp_monocots' => '#BBBB00',
		   '4_basals' => '#0000ff',
		   '4_C4_monocots' => '#000000', 
		   '7_C3_monocots' => '#FF0000',
		   '9_C3_dicots' => '#FF7700',
		  );
my %id_species = ();
my %id_color = ();



my %species_color = ();
while (my ($grp, $species) = each %$groups) {
#  print "$grp \n";
  for my $sp (keys %$species) {

    my $color = $group_color{$grp};
 #   print "$sp   $color \n";
    $species_color{$sp} = "&!color=$color";
    #   print "$sp       $grp \n";
  }
}
#exit;
# if (defined $gg_filename and -f $gg_filename) {
#   my $id_species = store_gg_info($gg_filename)
# } else {
#   warn "gg file: [$gg_filename] not defined or file doesnt exist. Use species in newick.\n";
# }

my @ids = ();

my $newick_expression = '';
if (defined $newicks_filename and -f $newicks_filename) {
  open my $fh_newick, "<", "$newicks_filename" or die "file [$newicks_filename] couldnt be opened for reading.\n";
  my @newick_lines = <$fh_newick>;
  while (my $line = shift @newick_lines) {
    # if ($line =~ /^Id M/) {
    #   #   print "$line  $id  \n";
    # }
    if ($line =~ /^Id (\S+)/ and $1 eq $id) {
      #  $id_color{$id} = 
      #     print "AAA: $line \n";
      my $nwck_line = shift @newick_lines;
      #      print "BBB: $nwck_line \n";
      if ($nwck_line =~ /^\s*\S+\s+(.*)/) {
	$newick_expression = $1;
	last;
      }
    }
  }
} else {
  die "newicks file: [$newicks_filename] no defined or file doesnt exist.\n";
}
$newick_expression =~ s/\s//g; # remove whitespace from newick expression.
#print "$newick_expression \n";

 # while(my($s, $c) = each %species_color){
 # print "$s $c \n";
 # }
 # exit;

#         [&!color=#02f6f0]
my $n_leaves = 0;
my $count = 0;
while (1) {
  if ($newick_expression =~ s/([,(])([^(][^[]+)[[]species=([^]]+)[]]/$1dummy_seqid[color=dummy_color]/) {
    $n_leaves++;
    my $seq_id = $2;
    my $species = $3;
    my $color = $species_color{$species};
    #   $id_species{$seq_id} = $species;

 #   print STDERR "$seq_id  $species  $color \n";
    my $seqid_w_sp = $seq_id . "__$species";
    if ($seq_id eq $id) {
      $color = "&!color=$qcolor"; # query color
    }
    $id_color{$seqid_w_sp} = $color;
 #   print STDERR "$seqid_w_sp  $color \n";
    #   if($id eq $seq_id) { print "$id  $seq_id $seqid_w_sp  $color  $qcolor  ", $id_color{$seqid_w_sp}, "\n"; }
    push @ids, $seqid_w_sp;
#       print "$seq_id  $seqid_w_sp   $species   $color \n";
    #    my $seqid_sp_string = $seq_id . "[$species]";
    $newick_expression =~ s/dummy_seqid/$seqid_w_sp/;
    #   $newick_expression =~ s/color=dummy_color/$color/;
    # $newick_expression =~ s/[[]color=dummy_color[]]//;
    $newick_expression =~ s/color=dummy_color/species = $species/;
    #print "$newick_expression \n";
    $count++;
    #exit if($count > 2);
  } else {
    last;
  }
}
#print "XXXXX \n";
#exit;
#### OUTPUT ####


my $space = '    ';
print "#NEXUS\n" , 
  "begin taxa;\n", 
  $space, "dimensions ntax=", $n_leaves, ";\n",
  $space, "taxlabels\n";
for (@ids) {
#print STDERR "Id: $_ \n";
  print $space, "'", $_, "'[", $id_color{$_}, "]\n";
}
print ";\n", "end;\n\n";

print "begin trees; \n",
  "$space", "tree my_tree = [&R] ";
$newick_expression =~ s/\][0-9.]+:([0-9.]+)/]:$1/g; # remove branch supports for terminal branches (due to rerooting I guess).
#print "XXXX $newick_expression;\n";
if ($keep_branch_supports) {
  $newick_expression =~ s/\)([0-9.]+):([0-9.]+)/)[branch_length=$1]:$2/g; # label the branch support values, i.e. like [branch_support=0.887]:0.345,
} else {
  $newick_expression =~ s/\)[0-9.]+:([0-9.]+)/):$1/g; # remove branch supports for other branches.
}
$newick_expression =~ s/\s+//g;
print "$newick_expression;\n";
#exit;
print "end;\n\n";

# open my $fh_ft, "<", "ftblock.txt" or die "XXXX\n";
# my @ftblock = <$fh_ft>;
# print join('', @ftblock), "\n";

# not clear this figtree block is needed:
# print $figtree_block;











sub store_gg_info {	 #xx    # read in gene-genome association file
  # store in hash (keys: seqids, values: species). e.g. ('ATxxxx' => Arabidopsis_thaliana');
  my $gg_filename   = shift;
  my %seqid_species = ();
  if ( defined $gg_filename ) {
    if ( -f $gg_filename ) {
      open my $fh_gg, "<", "$gg_filename";
      while (<$fh_gg>) {
	my @cols = split( " ", $_ );
	my $species = shift @cols;
	$species =~ s/:$//;	# remove final colon if present.
	for (@cols) {
	  if ( exists $seqid_species{$_} ) {
	    warn "key $_ already stored with species: ",
	      $seqid_species{$_}, "\n";
	  } else {
	    $seqid_species{$_} = $species;
	  }
	}
      }
      close $fh_gg;
    } else {	    # done storing gg_file info in hash %seqid_species
      die "$gg_filename: no such file.\n";
    }
  } else {
    die "gg filename undefined. \n";
  }
  return \%seqid_species;
}
