#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);
use Getopt::Long;
use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;
use TomfyMisc qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring ';

my $gg_filename = undef;
GetOptions(
	   'gg_file=s'           => \$gg_filename, #
	 );
my $id_species = undef;
if(defined $gg_filename  and  -f $gg_filename){
$id_species = store_gg_info($gg_filename);
print STDERR "done storing id-species association information.\n";
}

my @files;
if (scalar @ARGV > 1) {
  @files = @ARGV;
} elsif (scalar @ARGV == 1){
  my $file_of_filenames = shift; 
  open my $fh0, "<", "$file_of_filenames" or die "couldnt open [$file_of_filenames] for reading.\n";

  while (<$fh0>) {
 #   print;
    chomp;
    s/\s*(\S+)\s*/$1/;
    push @files, $_;
  }
}else{
  die "Must give filenames as arguments, or file of filenames as arg.\n";
}
my $n_sets = scalar @files;

# print join("\n", @files), "\n";
# exit;

# arguments are 2 or more filenames, of files with ids in first col.
# outputs list of ids, with string (e.g. 1 0 1 1 ) indicating which files it was present in,
# i.e. which venn diagram region it belongs in.
# and number of ids in each of these regions.

print "Input files: \n", join("\n", @files), "\n";
my %id_count = ();
my @file_id_hashes = ();
for my $the_file (@files){ # for each input file, store id:count key:value pairs
	my $file_id_hash = {};
	open my $fh, "<", "$the_file";
	while(<$fh>){
	  next if(/^\s*#/);
          next if(/^\s*$/);
          my ($id, $idcount) = (0, 0);
          if(/^\s*(\S+)\s+(\d+)/){
             ($id, $idcount) = ($1, $2);
          }elsif(/^\s*(\S+)\s/){
             ($id, $idcount) = ($1, 1);
          }
      #    print STDERR "idcount: $idcount \n";
			$file_id_hash->{$id} += $idcount; # stores ids in the file; one hash for each file
			$id_count{$id} += $idcount; # stores ids and counts (summed over all files).
	}
	push @file_id_hashes, $file_id_hash;
}
my %vennregion_count = ();
#my %qsp_idcount = (); # e.g.  3 2 0 1 1 1    would mean this id occurs 3 times in file 1, 2 x in file 2, etc.
		# 	'0 0 0 0 ' => 0, '0 0 0 1 ' => 0, '0 0 1 0 ' => 0, '0 0 1 1 ' => 0,
# 			'0 1 0 0 ' => 0, '0 1 0 1 ' => 0, '0 1 1 0 ' => 0, '0 1 1 1 ' => 0,
# 			'1 0 0 0 ' => 0, '1 0 0 1 ' => 0, '1 0 1 0 ' => 0, '1 0 1 1 ' => 0,
# 			'1 1 0 0 ' => 0, '1 1 0 1 ' => 0, '1 1 1 0 ' => 0, '1 1 1 1 ' => 0,
# );

for my $id (keys %id_count){
	my $count = $id_count{$id};
	my $id_category_string = '';
        my $qsp_idcount_string = ''; # e.g.  3 2 0 1 1 1    would mean this id occurs 3 times in file 1, 2 x in file 2, etc
	for (@file_id_hashes){
           my $idcount = (exists $_->{$id})? $_->{$id} : 0;
        #   print STDERR "IDCOUNT: $idcount \n";
		$id_category_string .= (exists $_->{$id})? '1 ' : '0 ';
           $qsp_idcount_string .= "$idcount ";
	}
	$vennregion_count{$id_category_string}++;
	my $n_ones = sum(split(" ", $id_category_string)); # number of files this id occurs in.
        my $n_fams_this_id = sum(split(" ", $qsp_idcount_string)); # number of occurences of this id, summed over files
        my $species = (defined $id_species->{$id})? $id_species->{$id} : '-';
	print "$id   $id_category_string  $n_ones    $qsp_idcount_string  $n_fams_this_id    $species \n";
}

my @nregion_histogram = ((0) x $n_sets);
my @skeys = #sort keys %vennregion_count;
  sort { sum(split(" ", $a)) <=> sum(split(" ", $b)) } keys %vennregion_count;
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
