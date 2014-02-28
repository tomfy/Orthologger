#!/usr/bin/perl -w
use strict;
use Getopt::Long;
print "# select_from_cladesout.pl  ", join(" ", @ARGV), "\n";
my $nesting_string = ''; #
my $min_paralogs = 0;
my $max_paralogs_string = '';
my $max_disallowed_string = '';
my $clades_required_string = '';

# Process long cl options
GetOptions(
    'nesting=s'      => \$nesting_string, # e.g. '1<2, 2<3' means clade 1 nested strictly inside clade 2, and 2 strictly inside 3.
'minparalogs=i' => \$min_paralogs, 
    'maxparalogs=s'      => \$max_paralogs_string, # e.g. '1:0, 2:0' means no other seqs of same species as query, in clade 1, clade 2.
    'maxdisallowed=s' => \$max_disallowed_string, # e.g. '1:0, 2:0, 3:0'  (would imply no limit in clade 4)
    'clades_required=s' => \$clades_required_string, # e.g. '1,2,3,4'
);


my @clades_required = split('', $clades_required_string);
my @the_types = ('NJ', 'ML', 'NJ_BS', 'ML_BS');
my ($prev_id, $id) = (undef, '');
my $type_acc_count = {};

reset_type_acc_count($type_acc_count, \@the_types);

my $first_line = <>;
$first_line =~ s/^\s*#\s*//; # remove initial # and any surrounding whitespace

my @cladespecs = split(/\s*:\s*/, $first_line);
print "#  ", join("  ", @cladespecs), "\n";
my $n_clades = scalar @cladespecs;
my $columns_per_clade = 3;
while (<>) {
#print;
$columns_per_clade = $1 if(/clade (\d) number[s]{0,1} are given/);
#print "cols per clade: $columns_per_clade \n";
#print "cols per clade: $columns_per_clade \n";
  next if(/^\s*#/);		# skip comment lines
  my @cols = split(" ", $_);
  $id = shift @cols;
  if (defined $prev_id  and  $id ne $prev_id) {
    print result_summary_string($prev_id, $type_acc_count, \@the_types);
    reset_type_acc_count($type_acc_count, \@the_types);
  }
my $type = shift @cols;
my @depths = ();
my @nparalogs =();
my @ndisallowed = ();
my @sizes = ();
  for(my $i=0; $i<$n_clades; $i++){
    $depths[$i] = @cols[$columns_per_clade*$i];
    $nparalogs[$i] = @cols[$columns_per_clade*$i + 1] - 1;
    $ndisallowed[$i] = @cols[$columns_per_clade*$i + 2];
    $sizes[$i] = ($columns_per_clade > 3)? @cols[$columns_per_clade*$i + 3] : -1;
# print "AAA:  ", $i, "  ", $depths[$i], "  ", $nparalogs[$i], "  ", $ndisallowed[$i], "  ", $sizes[$i], "\n";
}
#print join(", ", @depths), "  nesting_string: $nesting_string \n";
 my $nesting_OK = check_nesting(\@depths, $nesting_string);
my $paralogs_OK = check_paralogs(\@nparalogs, $max_paralogs_string, $min_paralogs);
my $disallowed_OK = check_disallowed(\@ndisallowed, $max_disallowed_string);
  my $all_OK = ($nesting_OK and $paralogs_OK and $disallowed_OK);


$type_acc_count->{$type}++ if($all_OK);
# print "OK  acp nested nparalogs nbad: [$OK]  [$all_clades_present]  [$nested]  [$n_paralogs]  [$n_bad_species] \n";
  # print "$id   " , $required_clades_present, "  ", $nested? 1 : 0, "  $n_paralogs  $n_bad_species  ", $OK? 'ACC' : 'REJ', "\n";
  $prev_id = $id;
my $outstring = $_;
chomp $outstring;

# print "$outstring  [$nesting_OK] [$paralogs_OK] [$disallowed_OK] {", $all_OK, "  ", ($nesting_OK and $paralogs_OK and $disallowed_OK), "}\n";

#print result_summary_string($id, $type_acc_count, \@the_types);

}
print result_summary_string($id, $type_acc_count, \@the_types);

sub result_summary_string{
  my $id = shift;
  my $type_acc_count = shift;
  my $the_types = shift;
  my $string = sprintf("%18s   ", $id);
  #  print join(";; ", @$the_types), "\n";
  #print join(";", keys %$type_acc_count), "\n";
  for (@$the_types) {
    #    print "$_ \n";
    $string .=			# $_ . ",  " . 
      sprintf("%4i  ", ((exists $type_acc_count->{$_})? $type_acc_count->{$_} : 0));;
  }
  $string .= "\n";
}

sub reset_type_acc_count{
  my $type_acc_count = shift;
  my $the_types = shift;
  for (@$the_types) {
    $type_acc_count->{$_} = 0;
  }
}

sub check_nesting{
  my $depths = shift;
  my $nesting_string = shift;
  my @nesting = split(/\s*,\s*/, $nesting_string);
#print "depths ", join("; ", @$depths), " ;;; $nesting_string \n";
  for(@nesting){
    my ($inner, $outer) = split(/\s*[<]\s*/, $_);
#print "$inner,  $outer \n";
    $inner--; $outer--;
    return 0 if($depths->[$inner] < 0 or $depths->[$outer] < 0); # one of the relevant clades not found.
    return 0 if($depths->[$inner] >= $depths->[$outer]); # one of the nesting conditions not met.
  }
 # all nesting conditions satisfied 
  return 1;
}

sub check_paralogs{
  my $nparalogs = shift;
  my $max_paralogs_string = shift;
my $min_paralogs = shift;
$min_paralogs = 0 if(!defined $min_paralogs);
  my @max_paralogs = split(/\s*,\s*/, $max_paralogs_string);
  for(@max_paralogs){
    my ($clade, $max_plogs) = split(/\s*[:]\s*/, $_);
    return 0 if($nparalogs->[$clade-1] > $max_plogs); # too many paralogs in this clade
    return 0 if($nparalogs->[$clade-1] < $min_paralogs);
  }
 # all paralogs conditions satisfied 
  return 1;
}
 sub check_disallowed{
  my $ndisallowed = shift;
  my $max_disallowed_string = shift;
  my @max_disallowed = split(/\s*,\s*/, $max_disallowed_string);
  for(@max_disallowed){
    my ($clade, $max_disallow) = split(/\s*[:]\s*/, $_);
#print "XXXXXXXXXXXXXXXXXXXXXXXXX [$clade][$max_disallow] \n";

    return 0 if($ndisallowed->[$clade-1] > $max_disallow); # too many paralogs in this clade
  }
 # all disallowed conditions satisfied 
  return 1;
}

