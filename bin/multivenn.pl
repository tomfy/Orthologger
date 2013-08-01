#!/usr/bin/perl -w
use strict;

my @files = @ARGV;

print "Input files: ", join(", ", @files), "\n";
my %id_count = ();
my @file_id_hashes = ();
for my $the_file (@files){
	my $file_id_hash = {};
	open my $fh, "<", "$the_file";
	while(<$fh>){
		if(/^\s*(\S+)\s/){
			$file_id_hash->{$1}++;
			$id_count{$1}++;
		}
	}
	push @file_id_hashes, $file_id_hash;
}
my %vennregion_count = ();
for my $id (keys %id_count){
	my $count = $id_count{$id};
	my $id_category_string = '';
	for (@file_id_hashes){
		$id_category_string .= (exists $_->{$id})? '1 ' : '0 ';
	}
	$vennregion_count{$id_category_string}++;
	print "$id  $id_category_string \n";
}

my @skeys = sort keys %vennregion_count;
for (@skeys){
print "$_  ", $vennregion_count{$_}, "\n";
}