#!/usr/bin/perl -w
use strict;

# read in a list of query ids
# pick those newicks out of newicks file
# generate a figtree nexus file for each one
# ( with groups color coded. )

my $qids_filename_or_id = shift;	# has ids in first column
my $newicks_filename = shift;
my $with_dotn = shift || 0;
my $groups = shift || 'AM';     # or C4
my $all_trees = ($qids_filename_or_id eq 'ALL')? 1 : 0;


my %qids = ();
if (open my $fh_qids, "<", "$qids_filename_or_id") { #  or die "Cant open $qids_filename_or_id for reading. \n";
   while (<$fh_qids>) {
      if (/^\s*(\S+)/) {
         my $theid = $1;
         $theid =~ s/\.\d+\s*$// if(! $with_dotn);
         $qids{$theid}++;
      }
   }
   close $fh_qids;
} else {            # interpret $qids_filename_or_id as a single sequence id
   $qids{$qids_filename_or_id}++;
}
for (keys %qids) {
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
      if ($all_trees or exists $qids{$nogenemodel_id} or 
           exists $qids{$id}) {
         #$qids{$nogenemodel_id} += 100;
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
      } else {
         #  not one of the ids of interest.
      }
   }
}
for (keys %qids) {
   print "$_  ", $qids{$_}, "\n";
}
