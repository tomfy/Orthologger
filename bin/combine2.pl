#!/usr/bin/perl -w
use strict;

# two files each with ids in first col.
# for ids in both files, print info from both files

my $f1 = shift;
my $f2 = shift;

my $id_info1 = store($f1);
my $id_info2 = store($f2);

for my $id (keys %$id_info1){
	if(exists $id_info2->{$id}){
		print "$id  ", $id_info1->{$id}, "  ", $id_info2->{$id}, "\n";
	}
}


sub store{
	my $file = shift;
	my %id_info = ();
	open my $fh, "<", "$file" or die "Couldnt open $file for reading.\n";
	while(<$fh>){
		my @cols = split(" ", $_);
		my $id = shift @cols;
		my $info = join(" ", @cols);
		$id_info{$id} = $info;
	}
	return \%id_info;
}
