#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec qw( splitpath );

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}

use lib $libdir;
use TomfyMisc qw 'store_gg_info   newick_genspid2idgensp  newick_genspid2id__gensp  newick_idgensp2id__gensp  newick_idgensp2genspid  read_in_group_species read_in_group_color';

# read in trees (newick format),
# from a set of newick files matching the pattern
# and write out nexus format files for 
# use by figtree, with colors specified
# which depend on species
# species are determined either from 1) a genome-gene association file if specified, or
# from information in the newick file itself, either 2) in the format sequenceid[species=Genus_species] or
# 3) in the format Genus_species_sequenceid

my $gg_filename = undef;
my $newick_filename_pattern = undef;
my $group_species_file = undef;
my $group_color_file = undef;
my $default_color = '#444444';

my $keep_branch_supports = 1; # default: keep branch supports in newick expression
GetOptions(
	   'gg_file=s'      => \$gg_filename, # not needed if species info present in newicks file.
	   'pattern=s' => \$newick_filename_pattern,
           'group_species_file=s' => \$group_species_file,
           'color_file=s' => \$group_color_file,
	  );

my $newick_filenames_str = `ls $newick_filename_pattern`;
my @newick_filenames = split(" ", $newick_filenames_str);

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


my $id_species = (defined $gg_filename)? store_gg_info($gg_filename) : {};

my %species_color = ();

my $species_group = read_in_group_species($group_species_file);
my $group_color = read_in_group_color($group_color_file);
while (my($sp, $grp) = each %$species_group) {
   my $color = $group_color->{$grp} // $default_color;
   $species_color{$sp} = $color;
}

for my $newickfile (@newick_filenames) {
   my @ids = ();
   open my $fh_newick, "<", "$newickfile" or die "couldn't open $newickfile for reading.\n";
   my $newick_expression = <$fh_newick>;
   while(<$fh_newick>){ $newick_expression = $_; }
   close $fh_newick;
   $newick_expression =~ s/\s//g; # remove whitespace from newick expression.

   my @x = ($newick_expression =~ /\[species/g);
   my $n_leaves = scalar @x;
   if ($n_leaves == 0) { # not [species= ]format, convert to id__gensp format.
      ($newick_expression, $n_leaves) = newick_genspid2id__gensp($newick_expression, $id_species, \@ids);
   } else { # format is [species= ], convert to id__gensp format.
      ($newick_expression, $n_leaves) = newick_idgensp2id__gensp($newick_expression, $id_species, \@ids);
   }

   my ($v, $path, $output_filename) = File::Spec->splitpath($newickfile);
   $output_filename =~ s/(txt|newick)\s*$//;
   $output_filename .= 'nexus';
   open my $FHout, ">", "$output_filename";

   my $space = '    ';
   my $nexus_string = "#NEXUS\n" .
     "begin taxa;\n" .
       $space . "dimensions ntax=" . $n_leaves . ";\n" .
         $space . "taxlabels\n";

   for my $id (@ids) {
      my $species = $id_species->{$id};
      my $xid = $id .  "__" . $species; # . "[species=$species]";
      $nexus_string .=  $space . "'" . $xid . "'[" . '&!color=' 
        . ($species_color{$id_species->{$id}} // $default_color)
          . "]\n";
   }
   $nexus_string .= ";\n" . "end;\n\n";

   $nexus_string .= "begin trees; \n" .
     "$space" . "tree my_tree = [&R] ";
   $newick_expression =~ s/\][0-9.]+:([0-9.]+)/]:$1/g; # remove branch supports for terminal branches (due to rerooting I guess).
   if ($keep_branch_supports) {
      # label the branch support values, i.e. like [branch_support=0.887]:0.345,
      $newick_expression =~ s/\)([0-9.]+):([0-9.]+)/)[branch_support=$1]:$2/g; 
   } else {
      $newick_expression =~ s/\)[0-9.]+:([0-9.]+)/):$1/g; # remove branch supports for other branches.
   }
   $newick_expression =~ s/\s+//g;
   $nexus_string .=  "$newick_expression;\n";
   $nexus_string =~ s/;;/;/;
   $nexus_string .=  "end;\n\n";

   print $FHout "$nexus_string \n";
   close $FHout;

   print "Done writing nexus file: $output_filename.\n";

   # not clear this figtree block is needed:
   # print $figtree_block;
}                               # end loop over newick files


# sub store_gg_info {	 #xx    # read in gene-genome association file
#    # store in hash (keys: seqids, values: species). e.g. ('ATxxxx' => Arabidopsis_thaliana');
#    my $gg_filename   = shift;
#    my %seqid_species = ();
#    if ( defined $gg_filename ) {
#       if ( -f $gg_filename ) {
#          open my $fh_gg, "<", "$gg_filename" or die "Couldn't open $gg_filename for reading.\n";
#          while (<$fh_gg>) {
#             my @cols = split( " ", $_ );
#             my $species = shift @cols;
#             $species =~ s/:$//; # remove final colon if present.
#             for (@cols) {
#                if ( exists $seqid_species{$_} ) {
#                   warn "key $_ already stored with species: ",
#                     $seqid_species{$_}, "\n";
#                } else {
#                   $seqid_species{$_} = $species;
#                }
#             }
#          }
#          close $fh_gg;
#       } else {      # done storing gg_file info in hash %seqid_species
#          die "$gg_filename: no such file.\n";
#       }
#    } else {
#       die "gg filename undefined. \n";
#    }
#    return \%seqid_species;
# }


# sub genspid2idgensp{ # change format of species and id in newick expression, e.g.:
#    # Arabidopsis_thaliana_AT1g123456.3  ->  AT1g123456.1[species=Arabidopsis_thaliana]
#    my $newick_string = shift;
#    my $id_genusspecies = shift; # hashref
#    my $ids = shift;             # array ref
#    my $leaf_count = 0;
#    while ( $newick_string =~ /([(,])([A-Z][a-zA-Z]*_[a-zA-Z]+)_([^\s:]+)/ ) {
#       my ($lparenorcomma, $genus_species, $id) = ($1, $2, $3);
#       $id_genusspecies->{$id} = $genus_species;
#       push @$ids, $id;
#       $id =~ s/_/X___X/g;
#       $newick_string =~ s/[(,][A-Z][a-zA-Z]*_[a-zA-Z]+_[^\s:]+/$lparenorcomma$id [species= $genus_species ]/;
#       $leaf_count++;
#    }
#    $newick_string =~ s/X___X/_/g; # put the underscores back in the ids
#    $newick_string =~ s/\s+//g;
#    return ($newick_string, $leaf_count);
# }


# sub genspid2id__gensp{ 
# #  Arabidopsis_thaliana_AT1g123456.3 -> AT1g123456.3__Arabidopsis_thaliana
#    my $newick_string = shift;
#    my $id_genusspecies = shift; # hashref
#    my $ids = shift;             # array ref
#    my $leaf_count = 0;
#    my $newick_copy = $newick_string;
#    while ( $newick_copy =~ s/([(,])([A-Z][a-zA-Z]*_[a-zA-Z]+)_([^\s:,)]+)/$1 $3 __ $2/) {
#       my ($lparenorcomma, $genus_species, $id) = ($1, $2, $3);
#       $id_genusspecies->{$id} = $genus_species;
#       push @$ids, $id;
#       $leaf_count++;
#    }
#    $newick_copy =~ s/\s+//g;
#    return ($newick_copy, $leaf_count);
# }

# sub idgensp2id__gensp{ # change format of species and id in newick expression, e.g.:
#    #  AT1g123456.1[species=Arabidopsis_thaliana] -> AT1g123456.3__Arabidopsis_thaliana
#    my $newick_string = shift;
#    my $id_genusspecies = shift; # hashref
#    my $ids = shift;             # array ref
#    my $leaf_count = 0;
#    while ( $newick_string =~ s/([(,])([^[:,(\s]+)\[species=([^]]+)\]/$1 $2 __ $3/) {
#       my ($lparenorcomma, $id, $genus_species) = ($1, $2, $3); 
#       $id_genusspecies->{$id} = $genus_species;
#       push @$ids, $id;
#       $leaf_count++;
#    }
#    $newick_string =~ s/\s+//g;
#    return ($newick_string, $leaf_count);
# }


# sub idgensp2genspid{ # change format of species and id in newick expression, e.g.:
#    #  AT1g123456.1[species=Arabidopsis_thaliana] -> Arabidopsis_thaliana_AT1g123456.3
#    my $newick_string = shift;
#    my $id_genusspecies = shift; # hashref
#    my $ids = shift;             # array ref
#    my $leaf_count = 0;
#    while ( $newick_string =~ s/([(,])([^[]+)\[species=([^]]+)\]/$1$3_$2/) {
#       my ($lparenorcomma, $genus_species, $id) = ($1, $2, $3);
#       $id_genusspecies->{$id} = $genus_species;
#       push @$ids, $id;
#       $leaf_count++;
#    }
#    return ($newick_string, $leaf_count);
# }

# sub read_in_group_species{
#    my $group_species_filename = shift;
#    my %species_group = ();
#    if (defined $group_species_filename and -f $group_species_filename) {
#       open my $fhin, "<", $group_species_filename or die "Couldn't open $group_species_filename for reading.\n";
#       my $state = 0;
#       my $group = 'unknown';
#       my $species_in_group = '';

#       while (<$fhin>) {
#          $_ =~ s/[#].*$//;      # remove comments
#          if ($state == 0) {
#             if (s/^\s*\[//) {   # start of list of species
#                $state = 1;
#                $species_in_group .= $_; # anything after [ on the line goes into species_in_group
#             } elsif (/^\s*(\S+)/) {
#                $group = $1;
#             }
#          } elsif (/^\s*\]/) {
#             $state = 0;
#             $species_in_group =~ s/^\s+//;
#             $species_in_group =~ s/\s+$//;
#             my @species = split(/[,\s]+/, $species_in_group);
#             for my $sp (@species) {
#                $species_group{$sp} = $group;
#             }
#             $species_in_group = '';
#             $group = 'unknown';

#          } else {
#             $species_in_group .= $_;
#          }
#       }
#       close $fhin;
#    }
#    return \%species_group;
# }


# sub read_in_group_color{
#    my $group_color_filename = shift;
#    my %group_color = ();
#    if (defined $group_color_filename and -f $group_color_filename) {
#       open my $fhin, "<", $group_color_filename or die "Couldn't open $group_color_filename for reading.\n";
#       while (<$fhin>) {
#          next if(/^\s*#/);
#          /^\s*(\S+)\s+(\S+)/;
#          $group_color{$1} = $2;
#       }
#       close $fhin;
#    }
#    return \%group_color;
# }
