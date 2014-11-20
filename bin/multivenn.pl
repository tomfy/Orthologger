#!/usr/bin/perl -w
use List::Util qw(min max sum);
use strict;

my @files = @ARGV;

my $n_sets = scalar @files; 

# arguments are 2 or more filenames, of files with ids in first col.
# outputs list of ids, with string (e.g. 1 0 1 1 ) indicating which files it was present in,
# i.e. which venn diagram region it belongs in.
# and number of ids in each of these regions.

print "Input files: ", join(", ", @files), "\n";
my %id_count = ();
my @file_id_hashes = ();
for my $the_file (@files){ # for each input file, store id:count key:value pairs
	my $file_id_hash = {};
	open my $fh, "<", "$the_file";
	while(<$fh>){
	  next if(/^\s*#/);
		if(/^\s*(\S+)\s/){
			$file_id_hash->{$1}++; # stores ids in the file; one hash for each file
			$id_count{$1}++; # stores ids and counts (summed over all files).
		}
	}
	push @file_id_hashes, $file_id_hash;
}
my %vennregion_count = ();
		# 	'0 0 0 0 ' => 0, '0 0 0 1 ' => 0, '0 0 1 0 ' => 0, '0 0 1 1 ' => 0,
# 			'0 1 0 0 ' => 0, '0 1 0 1 ' => 0, '0 1 1 0 ' => 0, '0 1 1 1 ' => 0,
# 			'1 0 0 0 ' => 0, '1 0 0 1 ' => 0, '1 0 1 0 ' => 0, '1 0 1 1 ' => 0,
# 			'1 1 0 0 ' => 0, '1 1 0 1 ' => 0, '1 1 1 0 ' => 0, '1 1 1 1 ' => 0,
# );

for my $id (keys %id_count){
	my $count = $id_count{$id};
	my $id_category_string = '';
	for (@file_id_hashes){
		$id_category_string .= (exists $_->{$id})? '1 ' : '0 ';
	}
	$vennregion_count{$id_category_string}++;
	print "$id  $id_category_string \n";
}

my @nregion_histogram = ((0) x $n_sets);
my @skeys = sort keys %vennregion_count;
for (@skeys) {

 my $count = $vennregion_count{$_};
  $_ =~ s/\s+$//;
  my $n_region = sum(split(" ", $_));
  $nregion_histogram[$n_region-1] += $count;
  print "$n_region   $_   $count \n";
}

for (1..$n_sets) {
  print $_, "  ", $nregion_histogram[$_ - 1], "\n";
}

print "union size: ", sum(@nregion_histogram), "\n";


for my $i (0..$n_sets-1){
  for my $j ($i+1 .. $n_sets-1){
my $x = $files[$i] . "  " . $files[$j];
    my $idvennout = `idvenn.pl $x`;
chomp $idvennout;
print "$i  $j   $idvennout \n";
}
}
