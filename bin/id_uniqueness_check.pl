#!/usr/bin/perl -w
use strict;

# read two files (alfastas, or newicks, etc.)
# look for query ids which occur in both files

my $f1 = shift || '';
my $f2 = shift || '';

open my $fh1, "<", "$f1";
#open my $fh2, "<", "$f2";

my %id_count = ();

while(<$fh1>){
	next if(! /^Id (Med\S+)/);
	my $id = $1;
	$id_count{$id}++;
	if($id_count{$id} > 2){
		print "Id appears > 2 time in file $f1.  Id: $id  count: ", $id_count{$id}, "\n";	
	}
}
# print "size of id_count: ", scalar keys %id_count, "\n";
my @count_histogram = ();

while (my ($id, $count) = each %id_count){
#for my $id ( keys %id_count){
#	my $count = $id_count{$id};
#	print "$id  $count \n";
	$count_histogram[$count]++;
}

for (my $i=0; $i< scalar @count_histogram; $i++){
#	print "$i    $count_histogram[$i] \n";
}

if ( -f $f2){
	open my $fh2, "<", "$f2";
	while(<$fh2>){
		next if(! /^Id (Med\S+)/);
		my $id = $1;
	#	$id_count{$id}++;
		if(exists $id_count{$id}){
			print "Id appears in both files $f1, and $f2.  Id: $id  count: ", $id_count{$id}, "\n";
		}
	}
}

