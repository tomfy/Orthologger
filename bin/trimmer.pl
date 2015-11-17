#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use Devel::Cycle;
use List::Util qw (min max sum);

no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
use lib $libdir;

use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';

# read in a file with alignments for each of many families
# The format is for each family:
# 1) 1 line of info about family as whole
# 2) fasta for all aligned sequences in family (this is omitted in cases of families which do not
#    meet certain criteria; at present by default if they don't include sequences from >= 3 monocot species
# 3) one blank line
# for example: (this just shows the sequence (with gaps) for the first two of the
# Id Medtr1g004990.1 family. fam_size: 201 Arabidopsis_lyrata,Arabidopsis_thaliana,Brachypodium_distachyon,Brassica_rapa,Capsella_rubella,Carica_papaya,Chlamydomonas_reinhardtii,Cucumis_sativus,Glycine_max,Medicago_truncatula,Oryza_sativa,Physcomitrella_patens,Populus_trichocarpa,Ricinus_communis,Selaginella_moellendorffii,Solanum_lycopersicum,Solanum_tuberosum,Sorghum_bicolor,Thellungiella_halophila,Vitis_vinifera,Zea_mays  1

my ($input_fastas_file1, $input_fastas_file2);
my ($gt, $ct) = (0.6, 0.4);
my $output_fastas_file = 'output_fastas_file';


GetOptions(
	   'in1=s' => \$input_fastas_file1,
           'in2=s' => \$input_fastas_file2,
	   'consistency=f' => \$ct,
           'nongap_fraction=f' => \$gt,
	   'output_file=s' => \$output_fastas_file,
	  );

open my $fh1_in, "<", "$input_fastas_file1" or die "Couldn't open $input_fastas_file1 for reading.\n";
open my $fh2_in, "<", "$input_fastas_file2" or die "Couldn't open $input_fastas_file2 for reading.\n";

open my $fh_out, ">", "$output_fastas_file" or die "Couldn't open $output_fastas_file for writing.\n";


my $PID = $$;
my $comparesetfilename = 'comparesetfilename_' . $PID;
open my $fh_cset, ">", "$comparesetfilename";
my $tmpfile1 = $PID;
my $tmpfile2 = $PID;
$tmpfile1 .= '_tmp1'; # . "$input_fastas_file1";
$tmpfile2 .= '_tmp2'; # . "$input_fastas_file2";
print "# input alignment files:  $input_fastas_file1  $input_fastas_file2 \n";
print $fh_cset "$tmpfile1\n", "$tmpfile2\n";
close $fh_cset;

my $output_fastas_string = '';
while (1) {

   my ($idline1, $fasta1) = get_next_fasta_fam($fh1_in);
   my ($idline2, $fasta2) = get_next_fasta_fam($fh2_in);

   last if(!defined $idline1  or  !defined $idline2);

   my @cols1 = split(" ", $idline1);
   my @cols2 = split(" ", $idline2);

   my ($id1, $id2) = ($cols1[1], $cols2[1]);
   die "ids not equal; id1: $id1  id2: $id2  \n" if($id1 ne $id2);

   #print "id1, id2: $id1  $id2  \n";
   open my $fh1_out, ">", "$tmpfile1";
   open my $fh2_out, ">", "$tmpfile2";

   print $fh1_out "$fasta1 \n";
   print $fh2_out "$fasta2 \n";

   close $fh1_out;
   close $fh2_out;

my $tmp_out_fasta_filename = $PID . "_tmp_out.fasta";
   my $trimal_out = `trimal -compareset $comparesetfilename -gt $gt -ct $ct -out $tmp_out_fasta_filename`;
   # print "$trimal_out \n";
   $trimal_out =~ /File Selected:\s+(\S+)/;
   print STDERR $1, "\n";
   my $fasta_out_string = `cat $tmp_out_fasta_filename`;
   my $trimalout_idline = ($1 eq "$tmpfile1")? $idline1 : $idline2;
   $output_fastas_string .= "$trimalout_idline" . "$fasta_out_string" . "\n";
   print $fh_out $output_fastas_string;
   $output_fastas_string = '';
}

sub get_next_fasta_fam{
   my $fh_in = shift;
   my $idline = '';
   my $fasta = '';
   while (<$fh_in>) {
      if ( /^Id /) {
         my @cols = split( " ", $_ );
         my ( $qid, $fam_size, $taxa) = @cols[ 1, 4, 5];
         $idline = $_;
         $fasta  = '';
      } elsif ( /^\s*$/ ) {     # blank line -> process the sequence }
         return ($idline, $fasta);
      } else {
         $fasta .= $_;
      }
   }
   return (undef, undef);
}
