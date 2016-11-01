#!/usr/bin/perl -w
use strict;
use Getopt::Long;
no lib '/home/tomfy/cxgn/cxgn-corelibs/lib';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {     # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib'; 
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
use lib $libdir;

# read in file with similarity info for families (iis format)
# for each id get the corresponding sequence (fasta)
# usage example:
#  fam_iis2fasta.pl  -iis_file Mv21_1000.iis -fasta_in ../new21species-pep.fasta
# output (fasta sequences for each family) would be file:  Mv21_fam_1000.fastas if not specified with -out option

my $iis_file                  = undef; # blast output in iis format
my $input_fasta_filename      = undef; # fasta for all sequences
my $output_filename = undef;

# Process long cl options
GetOptions(
           'iis_file=s'          => \$iis_file,
           'fasta_infile=s'      => \$input_fasta_filename,
           'output_filename=s' => \$output_filename,
          );

if (!defined $output_filename) { # if no output filename specified on CL, make a file name
   $output_filename = $iis_file; # from $iis_file by replacing .iis with .fastas
   $output_filename =~ s/[.]iis$//;
   $output_filename .= '.fastas';
}
print STDERR "fams_iis2fasta.pl OUTPUT FILENAME: $output_filename \n";
########

my $id_sequence_all = store_fasta_sequences($input_fasta_filename);

open my $fh_iis, "<", "$iis_file" or die "couldnt open $iis_file for reading. \n";
open my $fh_out, ">", "$output_filename" or die "Failed to open $output_filename for writing.\n";
my ( $fam_size, $fam_string_head, $fam_string_fasta ) = ( 0, '', '' );
my %taxon_count = ();
my %id_present = ();            #
while ( my $line = <$fh_iis> ) {
   next if($line =~ /^\s*#/);   # skip comment lines
   next if($line =~ /^\s*$/);   # skip all whitespace lines
   if ($line =~ /^(\S+)\s*$/) {
      my @cluster_qids = split(",", $1);
      print $fh_out $fam_string_fasta, "\n" if($fam_string_fasta ne '');
      print $fh_out "Ids $1\n";
      $fam_string_fasta = '';
      %id_present = ();
      for my $qid (@cluster_qids) {
         $id_present{$qid} = 1;
         $fam_string_fasta .= ">$qid\n" . $id_sequence_all->{$qid} . "\n";
      }
   } elsif ($line =~ /^\s+(\S+)\s+(\S+)\s*$/) {
      my ($id2, $sim) = ($1, $2);
      $fam_string_fasta .= ">$id2\n" . $id_sequence_all->{$id2} . "\n" unless(exists $id_present{$id2});
   } else {
      warn "line has unexpected format: $line \n";
   }
}                              # loop over lines of blast output (iis)
print $fh_out $fam_string_fasta, "\n";


sub store_fasta_sequences {
   my $fasta_filename = shift;
   my %id_sequence    = ();

   open my $fh, "<", "$fasta_filename" or die "In store_fasta_sequences. Couldnt open $fasta_filename for reading.\n";
   my ( $id, $sequence ) = ( undef, '' );
   while ( my $line = <$fh> ) {
      if ( $line =~ /^>(\S+)/ ) {
         if ( defined $id ) {
            $sequence =~ s/[*]\s*$//;  # remove asterisk at end of sequence, if present.
            $id_sequence{$id} = $sequence;
         }
         $id       = $1;
         $sequence = '';
      } else {
         $line =~ s/^\s+//;     # remove initial whitespace
         $line =~ s/\s+$//;     # remove final whitespace
         $sequence .= $line;
      }
   }
   if ( defined $id ) {
      $sequence =~ s/[*]\s*$//;  # remove asterisk at end of sequence, if present.
      $id_sequence{$id} = $sequence;
   }
   return \%id_sequence;
}
