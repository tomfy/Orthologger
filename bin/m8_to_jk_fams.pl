#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw( min max sum );
use sort 'stable'; # guarantees a stable sort is used (perl uses mergesort)
# m8 -> abc
# for each query, take (in addition to query) take the top j matches, plus
# the best k matches of each species from the remaining matches.
my $default_n_top = 50;
my $default_n_of_species = 2;
my $min_ev = 1e-250; # if e-val is 0, set it to this, so it is affected by penalty
my $max_eval = 1e-4;
my $ggfilename = undef;
my $m8_filename = undef;
my $abc_filename = '';
my $n_top = $default_n_top;
my $n_of_species = $default_n_of_species;

# Process long cl options
GetOptions(
	   'input_filename=s' => \$m8_filename,
	   'output_filename=s' => \$abc_filename,
           'ggfilename=s' => \$ggfilename,

           #   'fam_size_limit=i' => \$fam_size_limit, # default size limit is n_species * 4
           'n_of_species=i' => \$n_of_species, 
           'n_top=i' => \$n_top,
           'max_eval=f' => \$max_eval,
	  );

die "No input filename given. Exiting. \n" if(!defined $m8_filename);
open my $fh_m8_in, "<", "$m8_filename";

if (!defined $abc_filename or $abc_filename eq '') {
   $abc_filename = $m8_filename;
   $abc_filename =~ s/[._]m8//;	# remove final .m8 (or _m8) if present
   $abc_filename .= '_fams.abc';
}
print STDOUT "$abc_filename\n";
my %species_present = ();
my %id_species = ();
if (defined $ggfilename  and  -f $ggfilename) { # read in id-species info
   open my $fhgg, "<", "$ggfilename";
   while (<$fhgg>) {
      my @ids = split(" ", $_);
      my $species = shift @ids;
      $species_present{$species} = 1;
      $species =~ s/:$//;       # remove final colon
      for (@ids) {
         $id_species{$_} = $species;
      }				#
   }
}

my $n_species = scalar keys %species_present; # number of species in analysis (may be fewer in a family)

print STDERR 
  "# input file: $m8_filename \n",
  "# output file: $abc_filename \n",
  "# ggfilename: $ggfilename \n",
  "# max_eval : $max_eval \n",
  "# n_top: $n_top \n",
  "# n_of_species: $n_of_species \n";

my $query_abc_string = '';      # id1 id1 mineval
my $top_matches_abc_string = '';
my $extra_matches_abc_string = '';
#my %id2_ev = ();

my $init_old_id1 = 'xxxxxx xxxxxx';
my $old_id1 = $init_old_id1;
my $old_idpair = 'xxxxxxx xxxxxxx' ;
#my $the_line = '';
my $matches_so_far = 0;
my $first = 1;
my %matchspecies_count = ();
#my %id2locus_ = (); # key is id2 with final '.1' or '_P01', etc. removed, value is 1
open my $fh_abc_out, ">", "$abc_filename";
my ($id1, $id2, $evalue);
while (my $the_line = <$fh_m8_in>) {
   next if($the_line =~ /^\s*#/); # skip comment lines
   next if($the_line =~ /^\s*$/); # skip blank lines
   my @cols = split(" ", $the_line);
   if (scalar @cols == 12) {	# input is consistent with m8 format
      ($id1, $id2, $evalue) = @cols[0,1,10];
   } else {
      die "File $m8_filename does not appear to have correct format (m8). Exiting.\n";
   }

   if ($id1 ne $old_id1 and !$first) { # this is first line of new family. So first, output the old family:

      print $fh_abc_out
        "$old_id1  $old_id1  $min_ev \n",
          "$top_matches_abc_string",
            "$extra_matches_abc_string";

      %matchspecies_count = ();
      $top_matches_abc_string = '';
      $extra_matches_abc_string = '';
      $matches_so_far = 0;
   }
   $old_id1 = $id1;
   $first = 0;
   next if($evalue > $max_eval); #  skip if e-value too big.
   my $idpair = "$id1 $id2";
   next if($idpair eq $old_idpair); #  skip lower quality match solutions for this id1 id2 pair.
   $old_idpair = $idpair;
   next if($id2 eq $id1);       # skip query - handled separately
   if ($matches_so_far < $n_top) {
      $top_matches_abc_string .= "$id1  $id2  $evalue \n";
      $matches_so_far++;
   } else {
      my $id2_species = $id_species{$id2};
      $matchspecies_count{$id2_species}++;
      if ($matchspecies_count{$id2_species} <= $n_of_species) {
         $extra_matches_abc_string .= "$id1  $id2  $evalue \n";
      }
   }
}                               # end loop over lines in .m8 file

  print $fh_abc_out
        "$id1  $id1  $min_ev \n",
          "$top_matches_abc_string",
            "$extra_matches_abc_string \n";


##################################################################################################
