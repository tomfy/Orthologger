#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);
use Getopt::Long;
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
use TomfyMisc qw 'store_gg_info timestring ';

my $fams_filename = undef;      # iis format
my $req_species_filename = undef;
my $n_sp_required = undef;
my $gg_filename = '/home/tomfy/Aug2015multispeciesquery/55set.gg';

GetOptions(
	   'fams_filename=s'           => \$fams_filename, #
	   'gg_filename=s'          => \$gg_filename,      # 
           'req_species_filename=s' => \$req_species_filename,
           'n_sp_required=i' => \$n_sp_required,
          );

# read in the set of species which we want to require some minimum
# number of -
my %reqsp = ();
open my $fh_sp, "<", "$req_species_filename" or die "couldn't open $req_species_filename for reading.\n";
while (<$fh_sp>) {
   next if(/^\s*#/);
   if (/^\s*(\S+)/) {
      $reqsp{$1} = 1;
   }
}

my $geneid_sp = store_gg_info("$gg_filename");

# read in the set of matches for each family, and output if there are sufficiently many
# species from the required set.
my %id2_sim = ();
my $out_string = '';
my $old_ids = 'not_an_id';
open my $fh_fams, "<", "$fams_filename" or die "couldn't open $fams_filename for reading.\n";
while (<$fh_fams>) {
   if (/^(\S+)/) {
      my $ids = $1;
      if ($old_ids ne 'not_an_id') {
         my %reqsp_present = ();
         my @qids = split(',', $old_ids);
         for my $qid (@qids) {
            my $sp = (exists $geneid_sp->{$qid})? $geneid_sp->{$qid} : undef;
            if (defined $sp) {
               $reqsp_present{$sp}++ if(exists $reqsp{$sp});
            } else {
               warn "species of seq $qid is unknown.\n";
            }
         }
         while (my($id2, $sim) = each %id2_sim) {
            my $sp = (exists $geneid_sp->{$id2})? $geneid_sp->{$id2} : undef;
            if(defined $sp){
            $reqsp_present{$sp}++ if(exists $reqsp{$sp});
         }else{
            warn "species of seq $id2 is unknown.\n";
         }
         }
         if (scalar keys %reqsp_present >= $n_sp_required) {
            print $out_string;
         }
      }
      %id2_sim = ();
      $out_string  = $_;
      $old_ids = $ids;
   } else {
      if (/^\s+(\S+)\s+(\S+)/) {
         $out_string .= $_;
         $id2_sim{$1} = $2;
      } else {
         warn "line has unexpected format: [[$_]] \n";
      }
   }
}
