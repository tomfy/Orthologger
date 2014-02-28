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
my $max_fam_size = 100;
my $ggfilename = undef;
my $first_of_species_factor = 1; # has no effect if no ggfile. # e20; #
# by setting this to some large number (e.g. 1e20) will keep the best match
# from each species, if it is within this factor of the e-value threshold
# so if it is 1, it has no effect.

# my $best_model_only = 0;

# Process long cl options
GetOptions(
	   'max_eval=s' => \$max_eval,
	   'fam_size_limit=i' => \$max_fam_size,
	   'ggfilename=s' => \$ggfilename
	   #	   'only_best_model=s' => \$best_model_only
	  );

my %id_species = ();
if(defined $ggfilename  and  -f $ggfilename){ # read in id-species info
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

# print STDERR "max_eval: $max_eval \n";

my $old_id1 = 'xxxxxx xxxxxx';
my $old_idpair = 'xxxxxxx xxxxxxx' ;
my $the_line = '';
my $id1_match_count = 0;
my $worst_eval_in_fam;
my %matchspecies_count = ();
#my %id2locus_ = (); # key is id2 with final '.1' or '_P01', etc. removed, value is 1
while ($the_line = <>) {
  my @cols = split(" ", $the_line);
  my $evcol = (scalar @cols == 3)? 2 : 10;
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
# print "$mxev $match_species $worst_eval_in_fam  $first_of_species_factor \n";
  #  my $id2trunc = $id2;
  #if(! $id2trunc =~ /scaffold
  #   $id2trunc =~ s/[.]\s+$//; # remove final .1  .2  etc. (arabidopsis and many others)
  # $id2trunc =~ s/_P\d{2}$//; # remove final _P01 etc. (maize)
  # $id2trunc = ~ s/[.]t\d+$//; # remove final .t1  etc. (beet)
  if ($evalue <= $mxev) {
    #	if(($idpair ne $old_idpair)  and  ($evalue <= $max_eval)){
    #$id2locus_infam{$id2trunc}++;
    print "$idpair $evalue \n"; # $worst_eval_in_fam ", scalar keys %matchspecies_count, "\n";	#$the_line;
    $matchspecies_count{$match_species}++;
    #		$old_idpair = $idpair;
  }
}
