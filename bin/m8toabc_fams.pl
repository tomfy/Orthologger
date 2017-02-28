#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw( min max sum );
use sort 'stable'; # guarantees a stable sort is used (perl uses mergesort)
# m8 -> abc
# a few different methods implemented ...
# 
my $default_max_of_a_species = 50;
my $default_n_top = 100;
my $min_ev = 1e-185; # if e-val is 0, set it to this, so it is affected by penalty
my $max_of_a_species = $default_max_of_a_species;
my $max_eval = 1e-4;
my $ggfilename = undef;
my $m8_filename = undef;
my $abc_filename = '';
my $fam_size_limit = undef;
my $penalty_type = 'quadratic'; # none, knee, quadratic
my $penalty_knee = 5; # beyond this number of matches of a species, penalize the e-value.
my $log10_eval_penalty = 10; #'inf'; # log10_eval_penalty = 'inf' => just take top n of each species 
my $diff1 = 7;               # 
my $diff2 = 2;       # $diff1 = 0, $diff2 = 4 -> log10(penalty) = 0,0,
my $species_list = [];          # array ref
my $n_top = $default_n_top;
my $max_eff_eval_factor = 1.01; # 1e200; # or set to just > 1

my %species_rank = (
# AM non-hosts:
                    'Arabidopsis_thaliana' => 1,
                    'Arabidopsis_lyrata' => 2,
                    'Capsella_rubella' => 3,
                    'Eutrema_salsugineum' => 4,
                    'Thellungiella_halophila' => 5,
                    'Brassica_rapa', => 6,
                    'Tarenaya_hassleriana' => 7,
                    'Spinacia_oleracea' => 8,
                    'Beta_vulgaris' => 9,
                    'Dianthus_caryophyllus' => 10,
                    'Utricularia_gibba' => 11,
                    'Nelumbo_nucifera' => 12,
                    'Spirodela_polyrhiza' => 13,
                    'Zostera_marina' => 14,

# non-angiosperms
                    'Volvox_carteri' => 1001,
                    'Ostreococcus_tauri' => 1002,
                    'Ostreococcus_lucimarinus' => 1003,
                    'Chlamydomonas_reinhardtii'  => 1004,
                    'Physcomitrella_patens'      => 1005,
                    'Selaginella_moellendorffii' => 1006,
                    'Picea_abies' => 1007,

# Amborella
                    'Amborella_trichopoda' => 108,

# AM host monocots
                    'Elaeis_guineensis' => 201, # oil palm (african)
                    'Phoenix_dactylifera'     => 202, # date palm
                    'Musa_acuminata' => 203,          # banana 
                    #                  'Setaria_italica'         => 1, # foxtail millet
                    'Triticum_urartu'       => 204, # wheat ancestor, diploid
                    'Phyllostachys_heterocycla' => 205, # bamboo
                    'Hordeum_vulgare'         => 206,   # barley
                    'Zea_mays'                => 207,   # maize
                    'Brachypodium_distachyon' => 208,
                    'Panicum_virgatum' => 209, # switchgrass
                    'Panicum_hallii' => 210,   #
                    'Sorghum_bicolor'         => 211,
                    'Oryza_sativa'            => 212, # rice

# AM host dicots:
                    'Aquilegia_coerulea' => 301, # columbine

                       'Solanum_lycopersicum' => 302, # tomato
                       'Solanum_tuberosum'    => 303, # potato
                       'Mimulus_guttatus' => 304,     # monkeyflower
                       'Fraxinus_excelsior' => 305,   # Ash
                       'Sesamum_indicum' => 306,

                       'Vitis_vinifera'       => 307, # grape

                       'Glycine_max'          => 308, # soy
                       'Phaseolus_vulgaris' => 309,
                     #  'Lupinus_angustifolius' => 1,
                       'Lotus_japonicus' => 310,
                       'Medicago_truncatula' => 311,

                       'Populus_trichocarpa'  => 312, # poplar
                       'Ricinus_communis'     => 313, # castor
                       'Cucumis_sativus'      => 314, # cucumber
                       'Manihot_esculenta' => 315,
                       'Salix_purpurea' => 316,

                       'Theobroma_cacao' => 317,
                       'Carica_papaya' => 318,
                       'Eucalyptus_grandis' => 319,
                       'Gossypium_raimondii' => 320,
                       'Citrus_clementina' => 321,
                       'Citrus_sinensis' => 322,

                    'Lupinus_angustifolius' => 380,
);

# Process long cl options
GetOptions(
	   'input_filename=s' => \$m8_filename,
	   'output_filename=s' => \$abc_filename,
           'ggfilename=s' => \$ggfilename,

	   'fam_size_limit=i' => \$fam_size_limit, # default size limit is n_species * 4
           'max_of_species=i' => \$max_of_a_species,
           'n_top=i' => \$n_top,
           'max_eval=f' => \$max_eval,
           'max_eff_eval_factor=f' => \$max_eff_eval_factor,

           'p_type=s' => \$penalty_type,
           'p_slope=s' => \$log10_eval_penalty,
           'p_knee=i' => \$penalty_knee,

           'p_diff1=f' => \$diff1,
           'p_diff2=f' => \$diff2,
	  );

my $log10_max_eff_eval = log($max_eval * $max_eff_eval_factor)/log(10.0);


my @pows = (0);        # these are the logs (base10) of the penalties.
if ($penalty_type eq 'quadratic') {
   my @d1s = ($diff1);
   for my $j (1..$max_of_a_species) {
      push  @d1s, $d1s[-1]+$diff2;
   }
   for my $j (1..$max_of_a_species) {
      $pows[$j] = $pows[$j-1]+$d1s[$j-1];
   }
} elsif ($penalty_type eq 'knee') {
   for my $j (1..$max_of_a_species) {
      $pows[$j] = ($j > $penalty_knee)? $log10_eval_penalty * ($j - $penalty_knee) : 0;
   }
} else {                        # no penalty
   for my $j (1..$max_of_a_species) {
      $pows[$j] = 0;
   }
}

my $pow_description = ($log10_eval_penalty eq 'inf')? "e-value only as tie-breaker \n" : "# log_10 eval penalty factors: " . join(", ", @pows[0..25]) . "\n";
print STDERR $pow_description;

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
if (!defined $fam_size_limit) {
   $fam_size_limit = 4 * $n_species;
}

print STDERR 
  "# input file: $m8_filename \n",
  "# output file: $abc_filename \n",
  "# ggfilename: $ggfilename \n",
  "# max_eval : $max_eval \n",
  "# fam size limit: $fam_size_limit \n",
  "# max of a species: $max_of_a_species \n", 
  "# n_top: $n_top \n";

  my %id2_ev = ();
my $init_old_id1 = 'xxxxxx xxxxxx';
my $old_id1 = $init_old_id1;
my $old_idpair = 'xxxxxxx xxxxxxx' ;
my $the_line = '';
my $fstring = ''; # for each id1, $fstring is concat of lines, each with "id2 eval \n";
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
         if ($log10_eval_penalty eq 'inf') {
            ($fam_string, $total_count, $species_present) = 
              get_fam_string_inf($old_id1, $fstring, \%id_species, \%species_rank, $fam_size_limit, $max_of_a_species, $n_top); #, $_24_species);
         } else {
            ($fam_string, $total_count, $species_present) = 
              get_fam_string($old_id1, $fstring, \%id_species, $fam_size_limit, \@pows);
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
} # end loop over lines in .m8 file

my ($fam_string, $top_count, $ff_count, $total_count, $species_top, $species_ff, $species_present); 

if ($log10_eval_penalty eq 'inf') {
   ($fam_string, $total_count, $species_present) = 
     get_fam_string_inf($old_id1, $fstring, \%id_species, \%species_rank, $fam_size_limit, $max_of_a_species, $n_top); #, $_24_species);
} else {
   ($fam_string, $total_count, $species_present) = 
     get_fam_string($old_id1, $fstring, \%id_species, $fam_size_limit, \@pows);
}
print $fh_abc_out "$fam_string \n";
print STDERR "$old_id1   $total_count  $species_present \n";

##################################################################################################

sub get_fam_string{
   my $id1 = shift;
   my $id2eval_string = shift;
   my $id_species = shift;      # keys ids, values species
   my $fam_size_limit = shift;
   #  my $max_n_of_species = shift;
   my $powers = shift;
   #  my $small = 1e-180;

   my $max_n_of_species = scalar @$powers - 1;
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
         my $log10eveff = $log10ev + $powers->[min($sp_count, $max_n_of_species)-1];
         $id2_log10eveff{$id2} = $log10eveff;
      } else {
         die "Id $id2 has no species defined. \n";
      }
   }

   my @sorted_id2s = 
     #   sort { $id2_eveff{$a} <=> $id2_eveff{$b} } keys %id2_eveff;
     sort { $id2_log10eveff{$a} <=> $id2_log10eveff{$b} } keys %id2_log10eveff;

   %species_count = ();
   my $qqev = (defined $id2_ev{$id1})? $id2_ev{$id1} : $min_ev;

   $out_fam_string .= "$id1  $id1  $qqev  " . sprintf("%7.2f", log($qqev)/log(10.0)) . "\n"; # make the query be a match to itself, and it is first match in family 
   # even if the e-value isn't that good.
   my $fam_size_count = 1;
   for (@sorted_id2s) {
      next if($_ eq $id1);      # don't put the query in again!
      my $ev = $id2_ev{$_};
      my $log10ev_eff = $id2_log10eveff{$_};
      last if($log10ev_eff > $log10_max_eff_eval);
      $out_fam_string .= sprintf("%24s  %24s  %8s  %7.2f\n", $id1, $_, $ev, $log10ev_eff);
      # "$id1  $_  $ev  $log10ev_eff \n";
      $fam_size_count++;
      my $species = $id_species->{$_};
      $species_count{$species}++;
      last if($fam_size_count >= $fam_size_limit);
   }
   chomp $out_fam_string;
   return ($out_fam_string, $fam_size_count, scalar keys %species_count);
}

sub get_fam_string_inf{ # version with 'infinite' penalty, i.e. just take best of each species, then 2nd best of each species, etc.
   my $id1 = shift;
   my $id2eval_string = shift;
   my $id_species = shift; # hashref keys: seq ids, values: species (e.g. 'Arabidopsis_thaliana')
   my $sp_rank = shift; # hashref representing ordering of species (keys: species, values: ranks)
   my $size_limit = shift; # family size limit
   my $max_of_species = shift || $default_max_of_a_species; # include at most this number of matches from each species
my $n_top = shift; #
  my %selid_ev = (); # 
   my %species_matches = (); # keys: species; values: arrayrefs of match ids, in order from small e-value to large.
   my %id2_ev = ();
   my @id2ev_lines = split("\n", $id2eval_string);
my $match_rank = 0;
   for (@id2ev_lines) {
      my ($id2, $ev) = ($1, $2) if(/^\s*(\S+)\s+(\S+)/);
      $ev = $min_ev if($ev < $min_ev); # impose min ev (so none are zero)
      $id2_ev{$id2} = $ev;
      my $species = $id_species->{$id2} || warn "no species for id: $id2!!!\n";
      if (exists $species_matches{$species}) {
         push @{$species_matches{$species}}, $id2;
      } else {
         $species_matches{$species} = [$id2];
      }
      if($match_rank < $n_top){
     #    print STDERR "ZZ  $id1  $id2  $match_rank  $n_top \n";
         next if($id2 eq $id1);
         $selid_ev{$id2} = $ev;
         $match_rank++;
      }
   }
   my @species_list = sort { $sp_rank->{$a} <=> $sp_rank->{$b} } keys %species_matches; # species occurring in this family, ordered according to $sp_rank.

   my $qqev = (defined $id2_ev{$id1})? $id2_ev{$id1} : $min_ev;
#$selid_ev{$id1} = $qqev;
   my $fam_size_count = 1;
   my $round_count = 0;
   my $out_fam_string = "$id1  $id1  $qqev \n";
   while ($fam_size_count < $size_limit and $round_count < $max_of_species) {
      my %ith_round_ev = ();	#
      my @ith_round_ids = ();
      for my $a_species (@species_list) {
         next if(! defined $species_matches{$a_species});
         my $matches = $species_matches{$a_species};
         if (scalar @$matches > 0) {
            my $match_id = shift @$matches; # get the best remaining match of this species
            my $match_ev = $id2_ev{$match_id};
            $ith_round_ev{$match_id} = $match_ev;
            push @ith_round_ids, $match_id;
         } else {               # no more matches of this species
            delete $species_matches{$a_species};
         }
      }

      my @sms = sort {$ith_round_ev{$a} <=> $ith_round_ev{$b}} @ith_round_ids; # sort the matches in this round
      # stable sort, so in case of equal e-values, order in @ith_round_ids is preserved (and it is dictated by @species_list). So ordering is reproducible.
      for my $match_id (@sms) {
         next if($match_id eq $id1);
         my $the_ev = $ith_round_ev{$match_id};
       #  $out_fam_string .= "XX $id1  $match_id  $the_ev \n";
         $selid_ev{$match_id} = $the_ev;
         $fam_size_count++;
         last if(scalar keys %selid_ev >= $size_limit); # ($fam_size_count >= $size_limit);
      }
      $round_count++;
      last if(scalar keys %species_matches == 0);
   }
   
   my @sms = sort {$selid_ev{$a} <=> $selid_ev{$b}} keys %selid_ev;
   for my $m_id (@sms){
      $out_fam_string .= "$id1  $m_id  ". $selid_ev{$m_id}. "\n";
   }
   chomp $out_fam_string;
   return ($out_fam_string, 1 + scalar keys %selid_ev, scalar @species_list);
}
