#!/usr/bin/perl -w
use strict;

my $nesting = shift || 'ltltlt';
my $pclade = shift;
my $xbadclade = shift;
my $clades_required_string = shift;
if(!defined $pclade){ $pclade = 2; }
if(!defined $xbadclade){ $xbadclade = 3; }
if(!defined $clades_required_string){ $clades_required_string = '111'; }
my @clades_required = split('', $clades_required_string);
my @the_types = ('NJ', 'ML', 'NJ_BS', 'ML_BS');
my ($prev_id, $id) = (undef, '');
my $type_acccount = {};

reset_type_acccount($type_acccount, \@the_types);
while (<>) {
  next if(/^\s*#/);		# skip comment lines
  my @cols = split(" ", $_);
  $id = shift @cols;
  if (defined $prev_id  and  $id ne $prev_id) {
    print result_summary_string($prev_id, $type_acccount, \@the_types);
    reset_type_acccount($type_acccount, \@the_types);
  }
  my $type = 'ML';
  if (scalar @cols eq 10) {
    $type = shift @cols;
  }

  # my @clades_required = (1,1,0);
  my $required_clades_present = 1;
  for (my $i = 0; 3*$i < scalar @cols; $i += 1) {
    if ($clades_required[$i]  and  $cols[3*$i] < 0) {
      $required_clades_present = 0;
      last;
    }
  }
  if ($required_clades_present){
	#iprint ":::: ", $_;
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
	}elsif($nesting eq 'ltdcdc'){  # No Selaginella clade nesting required
	$nested = $cols[0] < $cols[3];
} else {			# ??? just use the default: ltltlt
    warn "Nesting option: $nesting not implemented. Using ltltlt.";
    $nested = (
	       ($cols[0] < $cols[3]) and
	       ($cols[3] < $cols[6])
	      );
  }

  my $n_paralogs = ($pclade > 0)? $cols[3*$pclade - 2] - 1 : 0; # paralogs defined as other medicago in the clade specified by pclade. pclade = 0 -> nothing counts as a paralog.
#  print "$pclade  pclade > 0: ", ($pclade > 0)? '1': '0', ".  n paralogs: $n_paralogs \n";
# sleep(0.5);
  my $n_bad_species = $cols[3*$xbadclade - 1]; # disallowed species in the Selaginella clade.
  my $OK = ($required_clades_present  and  $nested  and  ($n_paralogs == 0)  and  ($n_bad_species == 0));
# print "$nested  $n_paralogs  $n_bad_species \n";  
$type_acccount->{$type}++ if($OK);
# print "OK  acp nested nparalogs nbad: [$OK]  [$all_clades_present]  [$nested]  [$n_paralogs]  [$n_bad_species] \n";
  # print "$id   " , $required_clades_present, "  ", $nested? 1 : 0, "  $n_paralogs  $n_bad_species  ", $OK? 'ACC' : 'REJ', "\n";
  $prev_id = $id;
}
print result_summary_string($id, $type_acccount, \@the_types);

sub result_summary_string{
  my $id = shift;
  my $type_acccount = shift;
  my $the_types = shift;
  my $string = sprintf("%18s   ", $id);
  #  print join(";; ", @$the_types), "\n";
  #print join(";", keys %$type_acccount), "\n";
  for (@$the_types) {
    #    print "$_ \n";
    $string .=			# $_ . ",  " . 
      sprintf("%4i  ", ((exists $type_acccount->{$_})? $type_acccount->{$_} : 0));;
  }
  $string .= "\n";
}

sub reset_type_acccount{
  my $type_acccount = shift;
  my $the_types = shift;
  for (@$the_types) {
    $type_acccount->{$_} = 0;
  }
}
