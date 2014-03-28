#!/usr/bin/perl -w
use strict;

# read in an abc format blast output file
# i.e. each line has  id1 id1 e-value 
# and an query id.
# sort the matches and output them (abc format with query in col 0)
# to stdout
# usage: sortblastmatches.pl  'Medtr1g0049400.1' < xxx_blastout.abc

my $query_id = shift;
#print "query: $query_id \n";
my %idpair_eval = ();
while (<>) {
  my ($id1, $id2, $e_value) = split(" ", $_);
  if ($id1 eq $query_id) {
  #  print "$id1 $id2\n";
    my $id_pair = "$id1 $id2";
    if (exists $idpair_eval{$id_pair}) {
      if ($e_value < $idpair_eval{$id_pair}) {
	$idpair_eval{$id_pair} = $e_value;
      }
    }else{
      $idpair_eval{$id_pair} = $e_value;
    }
  }elsif ($id2 eq $query_id) {
  #  print "$id1 $id2\n";
    my $id_pair = "$id2 $id1";
    if (exists $idpair_eval{$id_pair}) {
      if ($e_value < $idpair_eval{$id_pair}) {
	$idpair_eval{$id_pair} = $e_value;
      }
    }else{
      $idpair_eval{$id_pair} = $e_value;
    }
  }
}
my @skeys = sort {$idpair_eval{$a} <=> $idpair_eval{$b}} keys %idpair_eval;
for (@skeys) {
  print "$_ ", $idpair_eval{$_}, "\n";
}
