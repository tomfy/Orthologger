#!/usr/bin/perl -w
use strict;

# read in a list of query ids
# pick those newicks out of newicks file
# generate a figtree nexus file for each one
# ( with groups color coded. )

my $qids_filename = shift;	# has ids in first column
my $newicks_filename = shift;
my $groups = shift || 'AM'; # or C4

my %qids = ();
open my $fh_qids, "<", "$qids_filename"  or die "Cant open $qids_filename for reading. \n";
while (<$fh_qids>) {
  if (/^\s*(\S+)/) {
    $qids{$1}++;
  }
}
close $fh_qids;
for(keys %qids){
  print "# $_   ", $qids{$_}, "\n";
}
#exit;

print STDERR scalar keys %qids, "\n";
print STDERR "[", join("] [", keys %qids), "]\n";
#exit;
open my $fh_newicks, "<", "$newicks_filename" or die "Cant open $newicks_filename for reading. \n";
while (<$fh_newicks>) {
  if (/^Id (\S+)\s/) {
    my $id = $1;
    print STDERR "$id \n";
    my $nogenemodel_id = $id;
    $nogenemodel_id =~ s/[.]\d{1,2}\s*$//;
    print STDERR $nogenemodel_id, "\n";
    if(exists $qids{$nogenemodel_id} or exists $qids{$id}){
$qids{$nogenemodel_id} += 100;
    #print "ID $id \n";
    my $newick_line = <$fh_newicks>;
    # print $newick_line, "\n";
    $newick_line =~ s/^\s*\S+\s*\(/(/; #remove stuff before first left paren.
    #	print $newick_line, "\n";
    open my $fhout, ">", "tmp.newick";
    print $fhout "Id $id \n", "XX ", $newick_line, "\n";
    close $fhout;
    my $out_nexus_filename = $id . ".nexus";
    system "newicks2figtreenexus.pl -new tmp.newick -id $id -group $groups > $out_nexus_filename ";
 }else{
    print "Id $nogenemodel_id not present in file $qids_filename \n";
 }
}
}
for (keys %qids){
   print "$_  ", $qids{$_}, "\n";
}
