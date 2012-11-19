#!/usr/bin/perl -w
use strict;

my $file_head = shift;
my $run1_file = $file_head . '.run1.p';
my $run2_file = $file_head . '.run2.p';

#`gnuplot -e "set multiplot"`;
my $command1 = "plot \'$run1_file\' using 1:2, \'$run2_file\' using 1:2;";
my $command2 = ""; # "plot 'b' using 1:2;";
my $command = "$command1 $command2";
my $gnuplot_command = "gnuplot -persist -e \"$command\"";

#$gnuplot_command .= "plot  'b' using 1:2";
print "$gnuplot_command \n";

my $gnuplotout = `$gnuplot_command`;
#`gnuplot -e "unset multiplot"`;
