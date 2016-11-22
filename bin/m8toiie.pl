#!/usr/bin/perl -w
use strict;

my $patterns_string = shift; # e.g. 'Mtr*.m8 Cpap*.m8' 
my @patterns = split(" ", $patterns_string);

my @filenames = ();
for my $pattern (@patterns){
   push @filenames,  split(" ", `ls $pattern`);
}
# print STDERR join("; ", @filenames), "\n";
# exit;

for my $filename (@filenames){
   open my $fhin, "<", "$filename" or die "Couldn't open $filename for reading.\n";
   print iie_string($fhin);
}

############################

sub iie_string{
   my $fh = shift;
   my $iie_str = '';
   my $first_line = <$fh>;
   my @cols = split(" ", $first_line);
   my ($qid, $id2, $ev) = @cols[0, 1, 10];
   $iie_str .= "$qid\n";
   $iie_str .= "  $id2 $ev\n";
   my ($old_qid, $old_id2) = ($qid, $id2);

   while ( my $line = <$fh>) {
      @cols = split(" ", $line);
      ($qid, $id2, $ev) = @cols[0,1,10];
      if ($qid ne $old_qid) {
         $iie_str .= "$qid\n";
         $iie_str .= "  $id2 $ev\n";
         $old_qid = $qid;
         $old_id2 = $id2;
      } elsif ($id2 ne $old_id2) {
         $iie_str .= "  $id2 $ev\n";
         $old_id2 = $id2;
      }else{
      #   $iie_str .= "  xxxx\n";
      }
   }
   return $iie_str;
}
