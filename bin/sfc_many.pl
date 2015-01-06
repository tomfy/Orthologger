#!/usr/bin/perl -w
use strict;

# Give it several directories.
# In each, there should be a file 'all.cladesout'
# process each with sfc_brsupp.pl

my $MinN_min_support = shift || 0.501;
my $DinN_min_support = shift || 0.501;
my $DinM_min_support = shift || 0.501;
my $MinNA_min_support = shift || 0.501;

my %id_nregioncount = ();

my $the_filename = 'all.cladesout';
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
	    '200_inf/muscle24/newicks_55_20/cladesout/',
	    '200_10/muscle4/newicks/cladesout/',
	    '200_10/muscle4/newicks_55_20/cladesout/',
	   );
if (0) {
  $base_dir = '/home/tomfy/Genome_data/50species/compact/ge20/ok702/200_10/';
  @dirs = (
	   '/mafft-best/newicks/multiFT/cladesout/',
	   '/mafft-best/newicks_55_20/multiFT/cladesout/', 
	   '/muscle24/newicks/multiFT/cladesout/',
	   '/muscle24/newicks_55_20/multiFT/cladesout/'
	  );
}
if (0) {
  $base_dir = '/home/tomfy/Genome_data/50species/compact/ge20/ok702/200_10/';
  @dirs = (
	   '/mafft-best/newicks/cladesout/',
	   '/mafft-best/newicks_55_20/cladesout/', 
	   '/muscle24/newicks/cladesout/',
	   '/muscle24/newicks_55_20/cladesout/'
	  );
}

if (1) {
  $base_dir = '/home/tomfy/Genome_data/50species/compact/ge20/ok791/200_10/';
  @dirs = (
	   '/mafft-best/ftphyml_newicks/cladesout/',
	   '/mafft-best/ftphyml_55_newicks/cladesout/', 
	   '/muscle24/ftphyml_newicks/cladesout/',
	   '/muscle24/ftphyml_55_newicks/cladesout/'
	  );
  $the_filename = 'all_tlr.cladesout';
}

if (0) {
  $base_dir = '/home/tomfy/Genome_data/50species/compact/ge20/ok791/200_10/';
  @dirs = (
	   '/mafft-best/ftphyml_newicks/cladesout/',
	   '/mafft-best/ftphyml_55_newicks/cladesout/', 
	   '/muscle24/ftphyml_newicks/cladesout/',
	   '/muscle24/ftphyml_55_newicks/cladesout/'
	  );
  $the_filename = 'all_r.cladesout';
}

#my $x = 0.5;
#print "# x min support values: $x $x $x $x \n";
print "# min support values:  $MinN_min_support $DinN_min_support $DinM_min_support $MinNA_min_support \n";
for (@dirs) {
  my $the_file = $base_dir . $_ . $the_filename;
  my $out_filename = $base_dir . $_ . 'all_ok_' 
    . $MinN_min_support . "_" 
      .  $DinN_min_support . "_" 
	. $DinM_min_support  . "_" 
	  . $MinNA_min_support . '.cladesout';
#  print STDERR "# file: $the_file \n";
  my $sfc_brsupp_outstring = 
    `sfc_brsupp.pl $MinN_min_support $DinN_min_support $DinM_min_support $MinNA_min_support < $the_file`;
  open my $fh_out, ">", "$out_filename";
  print $fh_out $sfc_brsupp_outstring, "\n";
  my @ok_lines = split("\n", $sfc_brsupp_outstring);
  print STDERR "# $the_file  ", scalar @ok_lines, " ok fams.\n";
  for (@ok_lines) {
    if (/^\s*(Med\S+)/) {
      my $id = $1;
      $id_nregioncount{$id}++;
    }
  }
}

my %nregion_population = ();

my @skeys = sort { $id_nregioncount{$a} <=> $id_nregioncount{$b} } keys %id_nregioncount;
for (@skeys) {
  my $count = (exists $id_nregioncount{$_})? $id_nregioncount{$_} : 0;
  print "$_   $count \n";
  $nregion_population{$count}++;
}
# print histogram
my $total = 0;
for (1..scalar @dirs) {
  my $pop = $nregion_population{$_};
  print "# $_   $pop\n";
  $total += $pop;
}

