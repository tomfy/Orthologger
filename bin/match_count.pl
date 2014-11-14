#!/usr/bin/perl -w
use strict;

my $eval_column = shift || 2;	# for abc, use 10 for m8
my $verbose = shift || 0;
# count the number of matches (one-way) for each id
# i.e. count how many times each id occurs in first col.
my %id1_count = ();
my %id2_count = ();

my %id1id2_eval = ();
my $init_id = 'IMPOSSIBLE - or highly improbable - ID';
my $old_id = $init_id;
my $id_count;
my $default_e_val_threshold = 0.000001;
my $e_val_threshold = shift || $default_e_val_threshold;
my $id_low_e_val_count = 1;
my $max_eval = -10000;
my ($cume_id_count, $cume_good_eval_id_count) = (0, 0);
my $query_count = 0;
while (<>) {
  my @cols = split(" ", $_);
  my ($id1, $id2, $eval) = ($cols[0], $cols[1], $cols[$eval_column]);
  $max_eval = $eval if($eval > $max_eval);
  my $idpair = "$id1 $id2";
  #	$id1id2_eval{$idpair} = $eval;
  $id1_count{$id1}++;
  $id2_count{$id2}++;
  if ($id1 eq $old_id) {
    #	print "$id1  $eval  $id_count  $id_low_e_val_count \n";
    $id_count++;
    $id_low_e_val_count++ if($eval <= $e_val_threshold);
  } else {		# print out info on old query id and count it.
    if($old_id ne $init_id){
      $query_count++;
      $cume_id_count += $id_count; $cume_good_eval_id_count += $id_low_e_val_count;
      print "$old_id  $id_count  $id_low_e_val_count  $query_count  $cume_id_count  $cume_good_eval_id_count \n" if($verbose);
    }
    $id_count = 1; $id_low_e_val_count = 1; $old_id = $id1;
	
  }
  #	print "$idpair ", $id1id2_eval{$idpair}, "\n" if( ($id_count % 10000) == 0);
}

$cume_id_count += $id_count; $cume_good_eval_id_count += $id_low_e_val_count;
print "$old_id  $id_count  $id_low_e_val_count \n" if($verbose); # ($old_id ne $init_id) and $verbose);
$query_count++;
print "# Query count: $query_count;  Max e-value: $max_eval  $cume_id_count  $cume_good_eval_id_count \n";
exit; 
for (keys %id1_count) {
  print STDERR "$_ ", $id1_count{$_}, "  ", $id2_count{$_}, "\n";
}

exit;
print scalar keys %id1id2_eval, "\n";
#exit;
for (keys %id1id2_eval) {
  my ($idone, $idtwo) = split(" ", $_);
  my $id21pair = "$idtwo $idone";
  my $eval21 = (exists $id1id2_eval{$id21pair})? $id1id2_eval{$id21pair} : '-';
  print "$_  ", $id1id2_eval{$_}, "  $eval21 \n";
}

