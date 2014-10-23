#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw( min max sum );
# m8 -> abc
# impose max e-value (i.e. keep only matches with e-value <= $max_eval)
# keep only best solution for each id1-id2 pair
# impose family size limit, but once reach fam size limit,
# also keep other matches as good as worst in fam
#

my $min_ev = 1e-180;

my $max_eval = 1e-4;
my $ggfilename = undef;
my $m8_filename = undef;
my $abc_filename = '';
my $fam_size_limit = 150;
my $fff = 'inf';

my $species_list = []; # array ref

# Process long cl options
GetOptions(
	   'input_filename=s' => \$m8_filename,
	   'output_filename=s' => \$abc_filename,
	   'fam_size_limit=i' => \$fam_size_limit,
	   'factor=f' => \$fff,
	   'ggfilename=s' => \$ggfilename,
	  );

my $max_of_a_species = 50;
my @factors = (1);
my @pows = (1);
for my $j (1..$max_of_a_species) {
  $pows[$j] = ($j-1)*$fff; # $fff * (0,1,2,3,4,  5,6,7,8,9, ... ) 'method A'
#  $pows[$j] = ($j <= 2)? 0: ($j-2)*($j-1) if($method eq 'AA'); # 0,0,2,6,12, 20,30,42,56,72,  ...
  $factors[$j] = ($pows[$j] < 2000)? 10**$pows[$j] : 1e2000;
}

print STDERR 
  "# input file: $m8_filename \n",
  "# output file: $abc_filename \n",
  "# ggfilename: $ggfilename \n",
  "# max_eval : $max_eval \n",
  "# fam size limit: $fam_size_limit \n",
  ($fff eq 'inf')? "e-value only as tie-breaker \n" : "# match handicap factors powers: " . join(", ", @pows[1..10]) . "\n";

die "No input filename given. Exiting. \n" if(!defined $m8_filename);
open my $fh_m8_in, "<", "$m8_filename";

if (!defined $abc_filename or $abc_filename eq '') {
  $abc_filename = $m8_filename;
  $abc_filename =~ s/[._]m8//;	# remove final .m8 (or _m8) if present
  $abc_filename .= '_fams.abc';
}
print STDOUT "$abc_filename\n";
my %id_species = ();
if (defined $ggfilename  and  -f $ggfilename) { # read in id-species info
  open my $fhgg, "<", "$ggfilename";
  while (<$fhgg>) {
    my @ids = split(" ", $_);
    my $species = shift @ids;
    $species =~ s/:$//;		# remove final colon
    for (@ids) {
      $id_species{$_} = $species;
    }				#
  }
}

my %id2_ev = ();
my $init_old_id1 = 'xxxxxx xxxxxx';
my $old_id1 = $init_old_id1;
my $old_idpair = 'xxxxxxx xxxxxxx' ;
my $the_line = '';
my $fstring = '';
my $first = 1;
my %matchspecies_count = ();
#my %id2locus_ = (); # key is id2 with final '.1' or '_P01', etc. removed, value is 1
open my $fh_abc_out, ">", "$abc_filename";
my ($id1, $id2, $evalue);
while ($the_line = <$fh_m8_in>) {
  next if($the_line =~ /^\s*#/); # skip comment lines
  next if($the_line =~ /^\s*$/); # skip blank lines
  my @cols = split(" ", $the_line);
  if (scalar @cols == 12) {	# input is consistent with m8 format
    ($id1, $id2, $evalue) = @cols[0,1,10];
  } else {
    die "File $m8_filename does not appear to have correct format (m8). Exiting.\n";
  }

  if ($id1 ne $old_id1 and !$first) { # this is first line of new family. So first, output the old family:
    %matchspecies_count = ();
    if ($fstring ne '') {
      my ($fam_string, $total_count,  $species_present);
      if($fff eq 'inf'){
	($fam_string, $total_count, $species_present) = 
	  get_fam_string_inf($old_id1, $fstring, \%id_species, $fam_size_limit); #, $_24_species);
      }else{
	($fam_string, $total_count, $species_present) = 
	  get_fam_string($old_id1, $fstring, \%id_species, $fam_size_limit, \@factors, \@pows);
      }

      print $fh_abc_out "$fam_string \n";
      print STDERR "$old_id1   $total_count  $species_present \n";
    }
    %id2_ev = ();
    $fstring = '';
  }
  $old_id1 = $id1;
  $first = 0;
  next if($evalue > $max_eval); #  skip if e-value too big.
  my $idpair = "$id1 $id2";
  next if($idpair eq $old_idpair); #  skip lower quality match solutions for this id1 id2 pair.
  $old_idpair = $idpair;
  $id2_ev{$id2} = $evalue;
  $fstring .= "$id2  $evalue \n";
}

 my ($fam_string, $top_count, $ff_count, $total_count, $species_top, $species_ff, $species_present); 

  if($fff eq 'inf'){
	($fam_string, $total_count, $species_present) = 
	  get_fam_string_inf($old_id1, $fstring, \%id_species, $fam_size_limit); #, $_24_species);
      }else{
	($fam_string, $total_count, $species_present) = 
	  get_fam_string($old_id1, $fstring, \%id_species, $fam_size_limit, \@factors, \@pows);
      }
print $fh_abc_out "$fam_string \n";
      print STDERR "$old_id1   $total_count  $species_present \n";

##################################################################################################

sub get_fam_string{
  my $id1 = shift;
  my $id2eval_string = shift;
  my $id_species = shift;
  my $fam_size_limit = shift;
  my $factors = shift;
  my $powers = shift;
#  my $small = 1e-180;

  my $max_n_of_species = scalar @$factors -1;
  my $out_fam_string = '';
  my %id2_log10eveff = ();
  my %id2_ev = ();
  my %species_count = ();
  my @id2ev_lines = split("\n", $id2eval_string);
  for (@id2ev_lines) {
    my ($id2, $ev) = ($1, $2) if(/^\s*(\S+)\s+(\S+)/);
    $ev = $min_ev if($ev < $min_ev);
    $id2_ev{$id2} = $ev;

    my $species = $id_species->{$id2} || warn "no species for id: $id2!!!\n";
    if (defined $species) {
      $species_count{$species}++;
      my $sp_count = $species_count{$species};
      my $log10ev = log($ev)/log(10.0);
      my $log10eveff = $log10ev + $powers->[min($sp_count, $max_n_of_species)];
      $id2_log10eveff{$id2} = $log10eveff;
    } else {
      die "Id $id2 has no species defined. \n";
    }
  }

  my @sorted_id2s = 
    #   sort { $id2_eveff{$a} <=> $id2_eveff{$b} } keys %id2_eveff;
    sort { $id2_log10eveff{$a} <=> $id2_log10eveff{$b} } keys %id2_log10eveff;

  my $got_query = 0;
  my $fam_size_count = 0;
  %species_count = ();
  for (@sorted_id2s) {
    my $ev = $id2_ev{$_};
    $out_fam_string .=  "$id1  $_  $ev \n";
    $fam_size_count++;
    $got_query = 1 if($_ eq $id1);
    my $species = $id_species->{$_};
    $species_count{$species}++;
    last if($fam_size_count >= $fam_size_limit);
  }
  if (! $got_query) { # make sure the query is in the family
    my $qqev = (defined $id2_ev{$id1})? $id2_ev{$id1} : $min_ev;
    $out_fam_string .= "$id1  $id1  $qqev \n";
  }
  chomp $out_fam_string;
  return ($out_fam_string, $fam_size_count, scalar keys %species_count);
}

sub get_fam_string_inf{
  my $id1 = shift;
  my $id2eval_string = shift;
  my $id_species = shift;
  my $size_limit = shift;	# family size limit

  #  my $first_few_count = 0;
  my $total_count = 0;
  #  my $top_count = 0;
  my $species_present_top = 0;
  my $species_present_ff = 0;
  my %species_matches = (); # keys: species; values: arrayrefs of match ids
  my %idpair_ev = ();
  my %id2_ev = ();
#  my $small = 1e-180;

  my @id2ev_lines = split("\n", $id2eval_string);
  for (@id2ev_lines) {
    my ($id2, $ev) = ($1, $2) if(/^\s*(\S+)\s+(\S+)/);
    $ev = $min_ev if($ev < $min_ev); # just for consistency with get_fam_string_A
    $id2_ev{$id2} = $ev;
    my $species = $id_species->{$id2} || warn "no species for id: $id2!!!\n";
    if (exists $species_matches{$species}) {
      push @{$species_matches{$species}}, $id2;
    } else {
      $species_matches{$species} = [$id2];
    }
  }
  my @species_list = keys %species_matches;

  my $fam_size_count = 0;
  my $got_query = 0;
  while ($fam_size_count < $size_limit) {
    my %ith_round_ev = ();	#
    for my $a_species (@species_list) {
      next if(! defined $species_matches{$a_species});
      my $matches = $species_matches{$a_species};
      if (scalar @$matches > 0) {
	my $match_id = shift @$matches;
	my $match_ev = $id2_ev{$match_id};
	$ith_round_ev{$match_id} = $match_ev;
      } else {			# no more matches of this species
	delete $species_matches{$a_species};
      }
    }

    my @sms = sort {$ith_round_ev{$a} <=> $ith_round_ev{$b}} keys %ith_round_ev;
    for my $match_id (@sms) {
      my $the_ev = $ith_round_ev{$match_id};
      $got_query = 1 if($match_id eq $id1);
      $idpair_ev{"$id1  $match_id"} = $id2_ev{$match_id};
      $fam_size_count++;
      last if($fam_size_count >= $size_limit);
    }
    last if(scalar keys %species_matches == 0);
  }
  if (! $got_query) {
    my $qqev = (defined $id2_ev{$id1})? $id2_ev{$id1} : $min_ev;
    $idpair_ev{"$id1  $id1"} = $qqev;
    $fam_size_count++;
  }

  my $out_fam_string = '';
  my @skeys = sort {$idpair_ev{$a} <=> $idpair_ev{$b}} keys %idpair_ev; # sort by e-value
  for (@skeys) {
    $out_fam_string .= "$_  " . $idpair_ev{$_} . "\n";
  }
  chomp $out_fam_string;
  return ($out_fam_string, $fam_size_count, scalar @species_list);
}


sub get_fam_string_old{
  # ($id1, $fstring, \%id_species, $fam_size_limit, $n_first_few, $ff_factor);
  my $id1 = shift;
  my $id2eval_string = shift;
  my $id_species = shift;
  my $fam_size_limit = shift;
  my $n_first_few = shift;
  my $ff_factor = shift;

  my $out_fam_string = '';
  my $first_few_count = 0;
  my $total_count = 0;
  my $top_count = 0;
  my $species_present_top = 0;
  my $species_present_ff = 0;
  my %species_count = ();
  # print STDERR "in get_fam_string. id2eval_string: \n $id2eval_string \n"; #exit;
  my @id2ev_lines = split("\n", $id2eval_string);
  my $top_n_max_eval; 
  for (my $i = 0; $i < min($fam_size_limit, scalar @id2ev_lines); $i++ ) {
    my $id2_ev = $id2ev_lines[$i];
    my ($id2, $ev) = ($1, $2) if($id2_ev =~ /^\s*(\S+)\s+(\S+)/);
    $top_n_max_eval = $ev;
    my $species = $id_species->{$id2} || undef;
    if (defined $species) {
      $species_count{$species}++;
      my $sp_count = $species_count{$species};
      $out_fam_string .= "$id1  $id2_ev \n";
      $total_count++;
      #     print STDERR "TOPN:  $id1  $id2   $species  $sp_count  $total_count \n";
    } else {
      die "Id $id2 has no species defined. \n";
    }
   
  }
  $species_present_top = scalar keys %species_count;
  $top_count = $total_count;
  if (scalar @id2ev_lines > $fam_size_limit) { # do any matches past the first $fam_size_limit
    for (my $i = $fam_size_limit; $i < scalar @id2ev_lines; $i++) {
      my $id2_ev = $id2ev_lines[$i];
      my ($id2, $ev) = ($1, $2) if($id2_ev =~ /^\s*(\S+)\s+(\S+)/);
      my $species = $id_species->{$id2} || undef;
      if (defined $species) {
	$species_count{$species}++;
	my $sp_count = $species_count{$species};
	if ( ( ($sp_count <= $n_first_few)  and  $ev <= ($top_n_max_eval * $ff_factor) ) or ($id1 eq $id2) ) {
	  $out_fam_string .= "$id1  $id2  $ev \n";
	  $first_few_count++;
	  #	  print STDERR "FF:    $id1  $id2   $species  $sp_count  $first_few_count \n";
	}
      } else {
	die "Id $id2 has no species defined. \n";
      }
    }
  }
  my $distinct_species_present = scalar keys %species_count;
  #print "[$out_fam_string]\n";
  $total_count += $first_few_count;
  chomp $out_fam_string;
  return ($out_fam_string, $top_count, $first_few_count, $total_count, $species_present_top, $distinct_species_present - $species_present_top, $distinct_species_present);
}


  # my $_24_species = [
  # 	       'Sorghum_bicolor', 'Zea_mays', 'Setaria_italica', 'Panicum_virgatum', # C4 grasses

  # 	       'Oryza_sativa', 'Hordeum_vulgare', 'Brachypodium_distachyon', 'Phyllostachys_heterocycla', # C3 grasses

  # 	       'Musa_acuminata', 'Phoenix_dactylifera', 'Spirodela_polyrhiza', # other C3 monocots

  # 	       'Amborella_trichopoda', # 4 basals
  # 	       'Picea_abies',
  # 	       'Selaginella_moellendorffii',
  # 	       'Physcomitrella_patens',

  # 	       'Aquilegia_coerulea', # columbine   # 9 C3 dicots
  # 	       'Solanum_lycopersicum',	    # tomato
  # 	       'Vitis_vinifera'      ,	    # grape
  # 	       'Medicago_truncatula',
  # 	       'Ricinus_communis'    ,	    # castor
  # 	       'Cucumis_sativus'     ,	    # cucumber
  # 	       'Arabidopsis_thaliana'   ,
  # 		   'Beta_vulgaris', # beet
  # 	       'Tarenaya_hassleriana',
  # 	      ];
