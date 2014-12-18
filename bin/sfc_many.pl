#!/usr/bin/perl -w
use strict;

# Give it several directories.
# In each, there should be a file 'all.cladesout'
# process each with sfc_brsupp.pl

my $MinN_min_support = shift || 0.25;
my $DinN_min_support = shift || 0.5;
my $DinM_min_support = shift || 0.4999;
my $MinNA_min_support = shift || 0.25;

my %id_nregioncount_x = ();
my %id_nregioncount_y = ();
my %id_nregioncount_xy = ();

my $base_dir = '/home/tomfy/Genome_data/50species/compact/ge20/';
my @dirs = (
	    '150_inf/mafftquick/newicks/cladesout/',
	    '200_inf/mafftquick/newicks/cladesout/',
	    '250_inf/mafftquick/newicks/cladesout/',

	    '150_inf/mafftquick/newicks_55_20/cladesout/',
	    '200_inf/mafftquick/newicks_55_20/cladesout/',
	    '250_inf/mafftquick/newicks_55_20/cladesout/',

	    '150_10/mafftquick/newicks/cladesout/',
	    '200_10/mafftquick/newicks/cladesout/',
	    '250_10/mafftquick/newicks/cladesout/',

	    '150_10/mafftquick/newicks_55_20/cladesout/',
	    '200_10/mafftquick/newicks_55_20/cladesout/',
	    '250_10/mafftquick/newicks_55_20/cladesout/',

	    '200_inf/muscle24/newicks/cladesout/',
	   );
my $x = 0.5;
for (@dirs) {
  my $the_file = $base_dir . $_ . 'all.cladesout';
  my $sfc_brsupp_outstring =  `sfc_brsupp.pl $x $x $x $x < $the_file`;
  my @ok_lines = split("\n", $sfc_brsupp_outstring);
  print STDERR "x: ", $the_file, "  ", scalar @ok_lines, " ok fams.\n";
  for (@ok_lines) {
    if (/^\s*(Med\S+)/) {
      my $id = $1;
      $id_nregioncount_x{$id}++;
      $id_nregioncount_xy{$id}++;
    }
  }

  $sfc_brsupp_outstring = `sfc_brsupp.pl $MinN_min_support $DinN_min_support $DinM_min_support $MinNA_min_support < $the_file`;
  @ok_lines = split("\n", $sfc_brsupp_outstring);
 print STDERR "y: ", $the_file, "  ", scalar @ok_lines, " ok fams.\n";
  for (@ok_lines) {
    if (/^\s*(Med\S+)/) {
      my $id = $1;
      $id_nregioncount_y{$id}++;
      $id_nregioncount_xy{$id}++;
    }
  }
}
#exit;
my @skeys = sort { $id_nregioncount_xy{$a} <=> $id_nregioncount_xy{$b} } keys %id_nregioncount_xy;
for (@skeys) {
  my $x_count = (exists $id_nregioncount_x{$_})? $id_nregioncount_x{$_} : 0;
 my $y_count = (exists $id_nregioncount_y{$_})? $id_nregioncount_y{$_} : 0;
  print "$_   $x_count  $y_count  " . $id_nregioncount_xy{$_} . "\n";
}
