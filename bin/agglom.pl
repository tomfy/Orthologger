#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw (min max sum);

my ($mu_bs_file, $mu_ph_file, $ma_bs_file, $ma_ph_file);

my %allid_count = ();
my %id_munjok = ();
my %id_muftok = ();
my %id_muphok = ();
my %id_munjbs = ();
my %id_muftbs = ();
my %id_manjok = ();
my %id_maftok = ();
my %id_maphok = ();
my %id_manjbs = ();
my %id_maftbs = ();

GetOptions(
	   'mu_bs_file=s' => \$mu_bs_file, # id  nj ft njbs40 ftbs40
	   'mu_ph_file=s' => \$mu_ph_file, # id nj ft ph 0 0 
	   'ma_bs_file=s' => \$ma_bs_file,
	   'ma_ph_file=s' => \$ma_ph_file,
	 );
my $fh;

open $fh, "<", "$mu_bs_file" or die "Could not open [$mu_bs_file] for reading.\n";
for (<$fh>) {
  if (/^\s*Med\S+/) {
    my @cols = split(" ", $_);

    my $id = $cols[0];
    $allid_count{$id}++;
    $id_munjok{$id} = $cols[1];
    $id_muftok{$id} = $cols[2];
    if (scalar @cols == 5) {
      $id_munjbs{$id} = $cols[3];
      $id_muftbs{$id} = $cols[4];
    } elsif (scalar @cols == 6) {
      $id_munjbs{$id} = $cols[4];
      $id_muftbs{$id} = $cols[5];
    } else {
      die "wrong number of columns: ", scalar @cols, "\n";
    }
  }
}
close $fh;

open $fh, "<", "$ma_bs_file" or die "Could not open [$ma_bs_file] for reading.\n";
for (<$fh>) {
  if (/^\s*Med\S+/) {
    my @cols = split(" ", $_);
    my $id = $cols[0];
    $allid_count{$id}++;
    $id_manjok{$id} = $cols[1];
    $id_maftok{$id} = $cols[2];
    if (scalar @cols == 5) {
      $id_manjbs{$id} = $cols[3];
      $id_maftbs{$id} = $cols[4];
    } elsif (scalar @cols == 6) {
      $id_manjbs{$id} = $cols[4];
      $id_maftbs{$id} = $cols[5];
    } else {
      die "wrong number of columns: ", scalar @cols, "\n";
    }
  }
}
close $fh;

open $fh, "<", "$mu_ph_file" or die "Could not open [$mu_ph_file] for reading.\n";
for (<$fh>) {
  if (/^\s*Med\S+/) {
    my @cols = split(" ", $_);
    my $id = $cols[0];
    $allid_count{$id}++;
    if ( (!exists $id_munjok{$id})  or  ($id_munjok{$id} != $cols[1]) ) {
      warn "Id: $id. Inconsistency in mu njok: ", $id_munjok{$id}, "  ", $cols[1];
      $id_munjok{$id} = $cols[1];
    }
    if (!exists $id_muftok{$id}  or  $id_muftok{$id} != $cols[2]) {
      warn "Id: $id. Inconsistency in mu ftok: ", $id_muftok{$id}, "  ", $cols[2];
    }
    $id_muphok{$id} = $cols[3];
  }
}
close $fh;

open $fh, "<", "$ma_ph_file" or die "Could not open [$ma_ph_file] for reading.\n";
for (<$fh>) {
  if (/^\s*Med\S+/) {
    my @cols = split(" ", $_);
    my $id = $cols[0];
    $allid_count{$id}++;
    if (!exists $id_manjok{$id}  or  $id_manjok{$id} != $cols[1]) {
      warn "Id: $id. Inconsistency in ma njok: ", $id_manjok{$id}, "  ", $cols[1];
      $id_manjok{$id} = $cols[1];
    }
    if (!exists $id_maftok{$id}  or  $id_maftok{$id} != $cols[2]) {
      warn "Id: $id. Inconsistency in ma ftok: ", $id_maftok{$id}, "  ", $cols[2];
    }
    $id_maphok{$id} = $cols[3];
  }
}
close $fh;

my $spacer = "  ";
for my $id (keys %allid_count){
$id_munjok{$id} = 0 if(!defined $id_munjok{$id});
$id_muftok{$id} = 0 if(!defined $id_muftok{$id});
$id_muphok{$id} = 0 if(!defined $id_muphok{$id});
$id_munjbs{$id} = 0 if(!defined $id_munjbs{$id});
$id_muftbs{$id} = 0 if(!defined $id_muftbs{$id});
$id_manjok{$id} = 0 if(!defined $id_manjok{$id});
$id_maftok{$id} = 0 if(!defined $id_maftok{$id});
$id_maphok{$id} = 0 if(!defined $id_maphok{$id});
$id_manjbs{$id} = 0 if(!defined $id_manjbs{$id});
$id_maftbs{$id} = 0 if(!defined $id_maftbs{$id});
}

for my $id (keys %allid_count) {
  my $mu_string = 
    $id_munjok{$id} . $spacer .
      $id_muftok{$id} . $spacer .
	$id_muphok{$id} . $spacer . $spacer .
	  $id_munjbs{$id} . $spacer .
	    $id_muftbs{$id} . $spacer;
  my $ma_string = 
    $id_manjok{$id} . $spacer .
      $id_maftok{$id} . $spacer .
	$id_maphok{$id} . $spacer . $spacer .
	  $id_manjbs{$id} . $spacer .
	    $id_maftbs{$id} . $spacer;

my $sumoks = $id_muftok{$id} + $id_muphok{$id} + $id_maftok{$id} + $id_maphok{$id};
my $sumbs = $id_muftbs{$id} + $id_maftbs{$id};

  print "$id   $mu_string     $ma_string    sumoks: $sumoks  sumbs: $sumbs \n", 
}
