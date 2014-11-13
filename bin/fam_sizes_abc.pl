#!/usr/bin/perl -w
use strict;

my %id_seqcount = ();

my $idpair_prev = 'xxx';
while(<>){
my @cols = split(" ", $_);

my $id1 = shift @cols;
my $id2 = shift @cols;
my $idpair = "$id1  $id2";
next if($idpair eq $idpair_prev);
$idpair_prev = $idpair;
# split(" ", $_);
$id_seqcount{$id1}++;
}

my $max_count = -1;
while (my ($id, $count) = each %id_seqcount){
	print "$id   $count \n";
	$max_count = $count if($count > $max_count);
} 
print "Max match count: $max_count \n";

