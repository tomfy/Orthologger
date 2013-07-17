#!/usr/bin/perl -w
use strict;

my $nesting = shift || 'ltltlt';
my $pclade = shift || 2;
my $xbadclade = shift || 3;

while (<>) {
  next if(/^\s*#/);		# skip comment lines
  my @cols = split(" ", $_);
  my $id = shift @cols;

my $all_clades_present = 1;
  for(my $i = 0; $i < scalar @cols; $i += 3){
    if($cols[$i] < 0){
      $all_clades_present = 0;
      last;
    }
  }
  my $nested;
  $nesting = 'ltltlt' if($nesting eq 'ltlt');
  if ($nesting eq 'ltltlt') {
    $nested = (
	       ($cols[0] < $cols[3]) and
	       ($cols[3] < $cols[6])
	      );
  } elsif ($nesting eq 'lelelt') {
    $nested = (
	       ($cols[0] <= $cols[3]) and
	       ($cols[3] <= $cols[6]) and
	       ($cols[0] < $cols[6])
	      );
  } elsif ($nesting eq 'lelele') {
    $nested = (
	       ($cols[0] <= $cols[3]) and
	       ($cols[3] <= $cols[6])
	      );
  } elsif ($nesting eq 'dcltlt') {
    $nested = (
	       ($cols[3] < $cols[6]) and
	       ($cols[0] < $cols[6])
	      );

  } elsif ($nesting eq 'dclelt') {
    $nested = (
	       ($cols[3] <= $cols[6]) and
	       ($cols[0] < $cols[6])
	      );

  } elsif ($nesting eq 'dclele') {
    $nested = (
	       ($cols[3] <= $cols[6]) and
	       ($cols[0] <= $cols[6])
	      );

  } elsif ($nesting eq 'dcdclt') {
    $nested = (
	       ($cols[0] < $cols[6])
	      );

  } elsif ($nesting eq 'dcdcle') {
    $nested = (
	       ($cols[0] <= $cols[6])
	      );
  } elsif ($nesting eq 'dcdcdc') { # No nesting requirement.
    $nested = 1;
  } else {			# ??? just use the default: ltltlt
    warn "Nesting option: $nesting not implemented. Using ltltlt.";
    $nested = (
	       ($cols[0] < $cols[3]) and
	       ($cols[3] < $cols[6])
	      );
  }

  my $n_paralogs = $cols[3*$pclade - 2] - 1; # paralogs defined as other medicago in the 3monocots clade. 
  my $n_bad_species = $cols[3*$xbadclade - 1]; # disallowed species in the Selaginella clade.
  my $OK = ($all_clades_present  and  $nested  and  ($n_paralogs == 0)  and  ($n_bad_species == 0));
  print "$id   " , $all_clades_present, "  ", $nested? 1 : 0, "  $n_paralogs  $n_bad_species  ", $OK? 'ACC' : 'REJ', "\n";
}
