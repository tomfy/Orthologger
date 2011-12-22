#!/usr/bin/perl -w
use strict;

# usage:  selclusters.pl selidfile < clusterfile
# output the clusters (same format as clusterfile)
# which contain a sequence whose id appears in selidfile

my $selid_file = shift;
open my $fhin, "<$selid_file";
my @selids = <$fhin>;
chomp foreach @selids;
my %selids = ();
foreach (@selids){
#print "selected id: [$_]\n";	
$selids{$_}++;
}

while(my $line = <>){
	my $clust_ids = $line;
	$clust_ids =~ s/^[^:]*[:]\s*//;
	my @ids = split(" ", $clust_ids);
	foreach my $id (@ids){
	$id =~ s/(.+)[(].*[)].*/$1/;		

#print $id, "\n"; 
if(exists $selids{$id}){
			print $line;
			last;
		}
	}
}
