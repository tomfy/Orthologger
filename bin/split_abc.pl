#!/usr/bin/perl -w
use strict;

# splits an abc file into pieces of approximately equal size.
# Does so in such a way that all the matches to one query is are in the
# same file, i.e. puts the breaks at points where the query id (leftmost column)
# changes.

my $input_abc_file = shift;
my $parts = shift || 4;
open my $fh_in, "<", "$input_abc_file";

my $part_number = 1;
my $prev_id1 = '-X_X_X_X_X_X_X-';
my $linesread = 0;
my $total_lines  = `wc $input_abc_file`;
$total_lines =~ /^\s*(\S+)/;
$total_lines = $1;
my $lines_per_part = int($total_lines/$parts + 1.0);

my $output_abc_file = part_filename($input_abc_file, $part_number);
#$output_abc_file =~ s/(part(\d+))([.]abc)/$1.$part_number$3/; # $input_abc_file . ".part$part_number";
#print "1: $1,  2: $2,  3: $3. \n";
my @part_filenames = ($output_abc_file);
open my $fh_out, ">", "$output_abc_file";
print STDERR "  $output_abc_file \n";


while (<$fh_in>) {
   #	print $fh_out $_;
   $linesread++;
   my @cols = split(" ", $_);
   my $id1 = $cols[0];
   if ($id1 ne $prev_id1) {
      if ($linesread >= $part_number * $lines_per_part) {
         close $fh_out;

         # print "$linesread  $output_abc_file \n";
         $part_number++;
         $output_abc_file = part_filename($input_abc_file, $part_number);
         print "$output_abc_file \n";
         push @part_filenames, $output_abc_file;
         #		$output_abc_file =~ s/(part(\d+))([.]abc)/$1.$part_number$3/; 
         open $fh_out, ">", "$output_abc_file";
         print STDERR "  $output_abc_file \n";
     
      }
      #	print $fh_out $_;
      $prev_id1 = $id1;
   }
   print $fh_out $_;
}
close $fh_out;
print STDOUT join(" ", @part_filenames), "\n";

sub part_filename{
my $filename = shift;
my $part_number = shift;
my $new_filename_end = '_part' . $part_number . '_fams.abc';
$filename =~ s/_fams[.]abc/$new_filename_end/;
return $filename;
}
