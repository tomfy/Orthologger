#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# m8 -> abc
# impose max e-value (i.e. keep only matches with e-value <= $max_eval)
# keep only best solution for each id1-id2 pair
# impose family size limit, but once reach fam size limit,
# also keep other matches as good as worst in fam
#

my $max_eval = 1e-12;		# default.
my $max_fam_size = 140;
my $ggfilename = undef;
my $m8_filename = undef;
my $abc_filename = undef;
my $first_of_species_factor = 1; # has no effect if no ggfile.
# by setting this to some large number (e.g. 1e20) will keep the best match
# from each species, if it is within this factor of the e-value threshold
# so if it is 1, it has no effect.

# my $best_model_only = 0;

# Process long cl options
GetOptions(
	   'input_filename=s' => \$m8_filename,
	   'output_filename=s' => \$abc_filename,
	   'max_eval=s' => \$max_eval,
	   'fam_size_limit=i' => \$max_fam_size,
	   'ggfilename=s' => \$ggfilename
	   #	   'only_best_model=s' => \$best_model_only
	  );

print STDERR "# input file: $m8_filename \n", "# output file: $abc_filename \n",
  "# max_eval: $max_eval \n", "# max_fam_size: $max_fam_size \n", "# ggfilename: $ggfilename \n";
die "No input filename given. Exiting. \n" if(!defined $m8_filename);
open my $fh_m8_in, "<", "$m8_filename";

if (!defined $abc_filename) {
  $abc_filename = $m8_filename;
  $abc_filename =~ s/[._]m8//;	# remove final .m8 (or _m8) if present
  $abc_filename .= '_fams.abc';
}

my %id_species = ();
if (defined $ggfilename  and  -f $ggfilename) { # read in id-species info
  open my $fhgg, "<", "$ggfilename";
  while (<$fhgg>) {
    my @ids = split(" ", $_);
    my $species = shift @ids;
    $species =~ s/:$//;		# remove final colon
    for (@ids) {
      $id_species{$_} = $species;
    }
  }
}

my $old_id1 = 'xxxxxx xxxxxx';
my $old_idpair = 'xxxxxxx xxxxxxx' ;
my $the_line = '';
my $id1_match_count = 0;
my $worst_eval_in_fam;
my %matchspecies_count = ();
#my %id2locus_ = (); # key is id2 with final '.1' or '_P01', etc. removed, value is 1
open my $fh_abc_out, ">", "$abc_filename";
while ($the_line = <$fh_m8_in>) {
  my @cols = split(" ", $the_line);
  my $evcol = undef;
  if (scalar @cols == 3) {
    $evcol = 2;
  } elsif (scalar @cols == 12) {
    $evcol = 10;
  } else {
    die "File $m8_filename does not appear to have correct format (m8 or abc). Exiting.\n";
  }
  my ($id1, $id2, $evalue) = @cols[0,1,$evcol];
  if ($id1 ne $old_id1) {
    $id1_match_count = 0;
    $old_id1 = $id1;
    $worst_eval_in_fam = $max_eval;
    %matchspecies_count = ();
  }
  next if($evalue > $max_eval); # skip if e-value too big.
  my $idpair = "$id1 $id2";
  next if($idpair eq $old_idpair); # skip lower quality match solutions
  $old_idpair = $idpair;
  $id1_match_count++;
  if ($id1_match_count == $max_fam_size) {
    $worst_eval_in_fam = $evalue;
  }
  my $match_species = (exists $id_species{$id2})? $id_species{$id2} : 'xxxxx xxxxx'; # if no ggfile, 'xxxxx xxxxx'
  my $mxev = (exists $matchspecies_count{$match_species})? $worst_eval_in_fam : $first_of_species_factor * $worst_eval_in_fam;
  if ($evalue <= $mxev) {
    # $id2locus_infam{$id2trunc}++;
    print $fh_abc_out "$idpair $evalue \n";
    $matchspecies_count{$match_species}++;
  }
}
