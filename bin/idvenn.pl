#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $file2 = shift;
# my $col = shift || 0;
my $verbose = shift || undef;
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
if($verbose){
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
}

my ($n1, $nboth, $n2) = (scalar @only1, scalar @both, scalar @only2);
print "# n 1, both, 2, union:  $n1   $nboth   $n2      ", $n1+$n2+$nboth, "\n"; 

#print "Union: \n\n", join("\n", (@only1, @both, @only2)), "\n";


sub store_ids{
	my $filename = shift;
	my %ids = ();
	open my $fh, "<", "$filename" or die "couldnt open file $filename \n";;
	while(<$fh>){
		next if(/^\s*#/); # skip commented lines
		$ids{$1}++ if(/^\s*(\S+)/); # store the id
	}
	return \%ids;
}
