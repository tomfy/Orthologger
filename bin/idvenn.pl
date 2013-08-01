#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $file2 = shift;

my $ids1 = store_ids($file1);
my $ids2 = store_ids($file2);

my @both = ();
my @only1 = ();
my @only2 = ();

for my $id1 (keys %$ids1){
	if(exists $ids2->{$id1}){
	push @both, $id1;
}else{
	push @only1, $id1;
}
}

for my $id2 (keys %$ids2){
	if(!exists $ids1->{$id2}){
	push @only2, $id2;
}
}

for(@only1){
  print $_, "  1only \n";
}
for(@both){
  print $_, "  both \n";
}
for(@only2){
  print $_, "  2only \n";
}
print "\n";

print "# n 1, both, 2: ", scalar @only1, "  ", scalar @both, "  ", scalar @only2, "\n"; 

#print "Union: \n\n", join("\n", (@only1, @both, @only2)), "\n";


sub store_ids{
	my $filename = shift;
	my %ids = ();
	open my $fh, "<", "$filename";
	while(<$fh>){
		next if(/^\s*#/);
		$ids{$1}++ if(/^\s*(\S+)/);
	}
	return \%ids;
}
