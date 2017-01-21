#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename 'dirname';
use File::Spec qw(splitpath);
use Cwd qw(abs_path getcwd);
use Time::Piece;

my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir); # collapses the bin/../lib to just lib
}
use lib $libdir;

use TomfyMisc qw(date_time);
use PhylogenomicPipeline;
use PipeSegGG;
use PipeSegBlast;


my %default_parameters = ( blast_max_matches => 2500, blast_max_e_value => 1e-6,
                           min_sequence_length => 20, n_pieces => 3,
                           blast_output_dir => 'blast_out',

                           families_output_dir => 'families',
                           family_size_limit => '500', family_max_e_value => 1e-6,
                           
                           family_fastas_output_dir => 'family_fastas',

                           alignment_output_dir => 'alignments',
                           alignment_program => 'mafft',
                           alignment_quality => 'quick',
                           trees_output_dir => 'newick_trees',
                           min_nongap_fraction => 0.7,
                           n_multi_ft => 1,
                           segments => {'init' => 1, 'blast' => 1, 'make_families' => 1, 'align' => 1},
                         );

my $control_filename = undef;
my $max_eval = 1e-6;		# default.
my $max_fam_size = undef;
my $n_pieces = 2;
my $blast_out_dir = 'blast_out';
my $fams_abc_dir = 'fams_abc';
my $fams_fastas_dir = 'fams_fastas';
my $mafft_alignment_dir = 'mafft_aligment_fastas';
my $muscle_alignment_dir = 'muscle_alignment_fastas';
my $n_multi_ft = 1;
#my $cl_query_taxon = undef;
my $query_id = undef;
my $default_min_sequence_length = 20;
my $query_number;               # either 'single' or 'multiple'
# my $added_groups_string = '';
my @pipeline_segments = ('gg', 'blast', 'make_families'); #, 'make_family_fastas', 'align_and_make_trees', 'end');
my %segment_class = (gg => 'PipeSegGG', blast => 'PipeSegBlast'); #, make_families => 'PipeSegMakeFams');
my %section_number = ();
my $i = 0;
for (@pipeline_segments) {
   $section_number{$_} = $i;
   $i++;
}
$section_number{'end'} = 100000;
my $start = undef;
my $last = undef;
my $default_start = 'blast';
my $default_last = 'end';

# Process long cl options
GetOptions(
	   'control_filename=s' => \$control_filename,

           'n_pieces=i' => \$n_pieces,
	   'max_eval=s' => \$max_eval,
	   'fam_size_limit=i' => \$max_fam_size,

           'n_multi_ft=i' => \$n_multi_ft,
           'query_id=s' => \$query_id,
  #         'query_taxon=s' => \$cl_query_taxon, # if specified on command line, overrides control file
           'start=s' => \$start, # 'blast', 'make_family', etc. 
           'last=s' => \$last,
	  );

die "Control filename must be specified on command line.\n" if(!defined $control_filename);
die "Specified control filename is: $control_filename; file does not exist.\n" if(! -f $control_filename);
my $pl_obj = PhylogenomicPipeline->new($control_filename, \%default_parameters);
print $pl_obj->stringify(), "\n\n";

# apply cl options overriding defaults & control file:
if (defined $start) {
   $pl_obj->set('start') = $start;
} else {
   $start = $pl_obj->get('start');
}

if (defined $last) {
   $pl_obj->set('last') = $last;
} else {
   $last = $pl_obj->get('last');
}


print STDERR "start and last segments: $start   $last \n";
die "problem with start or last: $start, $last \n" 
  if (!exists $section_number{$start} or !exists $section_number{$last} or ($section_number{$start}>$section_number{$last}));
#exit;

print STDERR "species, file: \n"; 
# exit;

# get date for incorporating into file names...
my $ltobj = localtime;
my $date = $ltobj->strftime('%Y%h%d');
my $time = $ltobj->hms;
#print STDERR "Pipeline start date, time: $date, ", $ltobj->hms, "\n\n";
print STDERR "Pipeline start date, time: $date, ", $ltobj->hms, "\n";
#exit;
my ($date, $time) = date_time();
print STDERR "Pipeline start date, time: $date, ", $ltobj->hms, "\n";
#exit;

my $qspecies = get_qspecies($pl_obj);
die "query files not specified.\n" if(!defined $qspecies);
# $qspecies = $species_long_short{$qspecies} if(exists $species_long_short{$qspecies});
print STDERR "qspecies: ", join(",", @$qspecies), "\n";
my $n_species = scalar keys %{$pl_obj->get('taxon_inputpath')};
my $filename_head = join('.', @$qspecies) . "_vs_" . "$n_species" . "_" . $date;
print STDERR "filename head: $filename_head \n";
my $progress_filename = $filename_head . '.progress';
print STDERR "progress_filename: $progress_filename \n";
open my $fh_progress, ">", "$progress_filename";
$pl_obj->set('fh_prog', $fh_progress);

my $gg_filename = $pl_obj->get('n_species') . "species_" . $date . ".gg";
my $all_species_fasta_filename = "all_" . $n_species . "_species_" . $date . ".fasta";
my @blast_out_m8_filenames = ();
my @fastas_filenames = ();


#my $blast_max_matches 
# $pl_obj->{blast_max_matches} = (exists $pl_obj->{blast_max_matches})? $pl_obj->{blast_max_matches} : 2500;
#my $blast_max_e_value = (exists $pl_obj->{blast_max_e_value})? $pl_obj->{blast_max_e_value} : '1e-8';

my $family_size_limit = (exists $pl_obj->{family_size_limit})? $pl_obj->{family_size_limit} : undef;
my $family_max_e_value = (exists $pl_obj->{family_max_e_value})? $pl_obj->{family_max_e_value} : '1e-12';
my $family_multiplicity_knee =  (exists $pl_obj->{family_multiplicity_knee})? $pl_obj->{family_multiplicity_knee} : 6;
my $family_log10_eval_penalty = (exists $pl_obj->{family_log10_eval_penalty})? $pl_obj->{family_log10_eval_penalty} : 12;
my $family_taxon_requirement = (exists $pl_obj->{family_taxon_requirement})? 
  $pl_obj->{family_taxon_requirement} : undef;
my $added_groups_string = (exists $pl_obj->{added_groups_string})?  $pl_obj->{added_groups_string} : undef;

my $alignment_program = (exists $pl_obj->{alignment_program})? $pl_obj->{alignment_program} : 'both'; # muscle, mafft, or both.
my $alignment_quality = (exists $pl_obj->{alignment_quality})? $pl_obj->{alignment_quality} : 'best'; # quick or best.

my $min_nongap_fraction = (exists $pl_obj->{min_nongap_fraction})? $pl_obj->{min_nongap_fraction} : 0.15;
my $min_seq_length = (defined $pl_obj->{min_sequence_length})? $pl_obj->{min_sequence_length} :  $default_min_sequence_length;

$n_multi_ft = $pl_obj->{n_multi_ft} if(defined $pl_obj->{n_multi_ft});
$n_multi_ft = ($n_multi_ft > 1)? $n_multi_ft : 1;
my $n_bs = ($n_multi_ft >= 1)? $n_multi_ft - 1 : 0;

print  $pl_obj->stringify(), "\n";
#exit;

  # ********** make gg string and file (gene-genome association)
my $Seg = 'PipeSegGG';
my $gg_obj;
#eval { $gg_obj = $Seg->new($pl_obj->get('taxon_inputpath'))};
my $evstr = '$gg_obj = ' . "$Seg" . '->new($pl_obj' . "->get('taxon_inputpath'))";
$evstr = '$gg_obj = ' . "$Seg" . '->new($pl_obj->get("taxon_inputpath"))';
eval '$gg_obj = ' . "$Seg" . '->new($pl_obj->get("taxon_inputpath"))';
print "evstr: $evstr \n";
#eval $evstr;
open my $fh_gg, ">", "$gg_filename" or die "couldn't open $gg_filename for writing.\n";
print $fh_gg $gg_obj->get('gg_string');
print $fh_progress "gg file: $gg_filename \n";
my $gg_summary = $gg_obj->summary();
print STDERR $gg_summary;
print $fh_progress $gg_summary;

#######################################################################################################
my $blast_obj = PipeSegBlast->new($pl_obj);
$blast_obj->run();
# die "exiting after blast.\n";
##########################################################################################


##########################################################################################
print STDERR $section_number{$start}, " ", $section_number{'make_families'}, " ", $section_number{$last}, "\n";
if (do_this_section('make_families', $start, $last, \%section_number)) {
   # *************** make family (abc) files *********************
   my @abc_part_filenames = ();
   print STDERR "blastout m8 filenames: ", join("; ", @blast_out_m8_filenames), "\n";
   for my $m8_filename (@blast_out_m8_filenames) {
      my ($v, $dir, $fname) = File::Spec->splitpath($m8_filename);
      $fname =~ s/m8$/abc/;
      my $fam_abcs_dir = $pl_obj->get('families_output_dir');
      mkdir $fam_abcs_dir unless(-d $fam_abcs_dir); # create dir if doesn't exist.
      my $fam_abcs_filename = $fam_abcs_dir . '/' . $fname;
      my $m8tofams_cl = "m8toabc_fams.pl -input $m8_filename -max_eval $family_max_e_value -gg $gg_filename -output $fam_abcs_filename ";
      $m8tofams_cl .= " -fam_size_limit $family_size_limit ";
      $m8tofams_cl .= "-multiplicity_knee $family_multiplicity_knee -log10penalty $family_log10_eval_penalty ";

      print STDERR "m8tofams_cl: $m8tofams_cl \n";
      my $abc_filename = `$m8tofams_cl`;
      # -fam_size_limit $family_size_limit `;
      $abc_filename =~ s/\s+$//; 
      print STDERR "abc_filename: [$abc_filename]\n";
      push @abc_part_filenames, $abc_filename;
   }
   #my $abc_part_filenames = `split_abc.pl $abc_filename $n_pieces`;
   #my @abc_parts = split(" ", $abc_part_filenames);
   #exit;

   print "[" ,join("][", @abc_part_filenames), "]\n";
   print $fh_progress "abc family filenames: \n", join("\n", @abc_part_filenames), "\n";

   # ************** done making family (abc) files *******************
die "exiting after making family (abc)files.\n";
   # ************* make family (fastas)

   die "fasta input file undefined or file does not exist.\n" unless(defined $all_species_fasta_filename and -f $all_species_fasta_filename);
   my $fam_fastas_dir = $pl_obj->{families_output_dir};
   mkdir $fam_fastas_dir unless(-d $fam_fastas_dir);
   for my $the_abc_part (@abc_part_filenames) {
      my ($v, $dir, $fname) = File::Spec->splitpath($the_abc_part);
      $fname =~ s/abc$/fastas/;
      my $output_fastas_filename = $fam_fastas_dir . '/' . $fname;
      #   print STDERR "ZZZZ: ", `which seq+matches2fasta.pl`, "\n";
      #  my $addgrpstring = "'" . 'amposdicots:Solanum_lycopersicum,Theobroma_cacao,Medicago_truncatula,Manihot_esculenta;monocots:Oryza_sativa,Phoenix_dactylifera' . "'";
      my $seqm2f_cl = "seq+matches2fasta.pl -gg $gg_filename -abc_file $the_abc_part -fasta_infile $all_species_fasta_filename -output_filename $output_fastas_filename";
      $seqm2f_cl .= " -taxon_requirement $family_taxon_requirement " if(defined $family_taxon_requirement);
      $seqm2f_cl .= " -added_groups $added_groups_string " if(defined $added_groups_string);
      print STDERR "seq+... cl:  $seqm2f_cl \n";
      system("$seqm2f_cl");
      # "seq+matches2fasta.pl -gg $gg_filename -abc_file $the_abc_part -fasta_infile $all_species_fasta_filename -output_filename $output_fastas_filename");
      push @fastas_filenames, $output_fastas_filename
   }

   print STDERR "Family fastas file ready to align: \n", join("\n", @fastas_filenames), "\n";

   print $fh_progress "Family_fastas_files: " , join("  ", @fastas_filenames), "\n";
}                               # end of make_families section
##########################################################################################
print "Quitting after make_families section.\n\n";
##########################################################################################
if (do_this_section('align_and_make_trees', $start, $last, \%section_number)) {
   # fork processes to do alignment, tree finding.
   my @alignment_programs = ('muscle', 'mafft'); # i.e. default is do both
   if ($alignment_program eq 'muscle') {
      @alignment_programs = ('muscle');
   } elsif ($alignment_program eq 'mafft') {
      @alignment_programs = ('mafft');
   }

   for my $align_program (@alignment_programs) {
      my $align_dir = $align_program . "_" . $pl_obj->{alignments_dir};
      mkdir $align_dir unless(-d $align_dir);
      die "Couldnt make directory $align_dir.\n" unless(-d $align_dir); # make a directory to put alignments in.
      chdir $align_dir or die "Couldnt change directory to $align_dir.\n";
      my $newick_dir = $pl_obj->{newick_trees_dir};
      mkdir $newick_dir unless(-d $newick_dir);
      my $n_alignment_files_to_do = scalar @fastas_filenames;
      my @alfastas_filenames = ();
      my @newicks_filenames = ();
      for my $a_fastas_filename (@fastas_filenames) {
         my $malign_out_filename = $a_fastas_filename;
         my ($v, $dir, $fname) = File::Spec->splitpath($malign_out_filename);
         $malign_out_filename = $fname;
         $malign_out_filename =~ s/fastas$/alfastas/;
         push @alfastas_filenames, $malign_out_filename;
         my $output_newick_filename = $newick_dir . "/" . $malign_out_filename;
         $output_newick_filename =~ s/alfastas$/newicks/;
         push @newicks_filenames, $output_newick_filename;
         my $pid = fork(); # returns 0 to child process, pid of child to parent process.
         if ($pid == 0) {  # child process
            $a_fastas_filename = '../' . $a_fastas_filename;
            $gg_filename = '../' . $gg_filename;
            print "just before malign. cwd: ", getcwd, "\n";
            my $malign_cl = "malign.pl  -input $a_fastas_filename  -align $align_program  -quality $alignment_quality  -output $malign_out_filename ";
            print STDERR "malign command line: [$malign_cl] \n";
            my $malign_stdout = `$malign_cl`;
            print "malign.pl finished aligning $a_fastas_filename; output file: $malign_out_filename. \n";

            my $tree_construct_cl =  "tree_construct.pl -gg $gg_filename -input $malign_out_filename -output $output_newick_filename -nongap_frac $min_nongap_fraction";
            # "NJFTPHtree.pl -gg $gg_filename -input $malign_out_filename -output $output_newick_filename -n_bs $n_bs -nongap_frac $min_nongap_fraction ";
            print "tree construction cl: $tree_construct_cl \n";
            my $tree_construct_stdout = `$tree_construct_cl`;
            exit(0);
         }
      }
      my $children_done = 0;
      while (wait() != -1) {
         $children_done++; print "Number of files finished aligning with $align_program: $children_done out of $n_alignment_files_to_do. \n";
      }

      my @files = split(" ", `ls *tmpfile`); # cleanup - delete temp files.
      for (@files) {
         unlink $_;
      }

      print "cwd: ", getcwd, "\n";
      chdir '../' or die "chdir ../ failed \n";
      print "cwd: ", getcwd, "\n";

      print $fh_progress "Done with aligning with $align_program and constructing trees. \n";
      print $fh_progress "alfastas files: \n" , join("\n", @alfastas_filenames), "\n";
      print $fh_progress "newick files: \n", join("\n", @newicks_filenames), "\n";
      print STDERR "Done with aligning with $align_program and constructing trees. \n";
      print STDERR "alfastas files: \n" , join("\n", @alfastas_filenames), "\n";
      print STDERR "newick files: \n", join("\n", @newicks_filenames), "\n";
   }                            # end of loop over alignment programs
}                               # end of align_and_make_trees section
close $fh_progress;
##########################################################################################

sub get_qspecies{
   my $plobj = shift;
   my $qsp_files = $plobj->get('query_taxon_inputpath');
  
   if ( (defined $qsp_files) and
        (ref $qsp_files eq 'HASH') ){
      my @sps = keys %{$qsp_files};
      return \@sps;
   }
   return undef;
}


# sub get_params_from_control_file{
#    my $control_filename = shift;
#    open my $fh_ctrl, "<", "$control_filename" or die "Could not open $control_filename for reading; exiting.";
#    my $added_groups_string = '';
#    my %param_name_val = ();
#    my %taxon_file = ();
#    while (<$fh_ctrl>) {
#       next if(/^\s*#/);         # skip all comment lines
#       if (/^\s*query_taxon_inputpath/) {
#          if (/^\s*query_taxon_inputpath\s+(\S+)\s+(\S+)/) { # $1 is taxon name, e.g. Zea_mays, $2 is path to corresponding input file (protein-models fasta format)
#             #    $taxon_file{$1} = $2;
#             $param_name_val{query_taxon} = $1;
#             $param_name_val{query_inputpath} = $2;
#          } elsif (/\S/) {       # not all whitespace
#             warn "A unexpected line in control_file: $_";
#          }
#       } elsif (/^\s*taxon_inputpath/) {
#          if (/^\s*taxon_inputpath\s+(\S+)\s+(\S+)/) {
#             my @exp_filenames = glob($2);
#             my $the_filename = $exp_filenames[0]; 
#             #   print STDERR "AAA: [", $2, "]  [", glob($2), "] [$the_filename]\n";            
#             $taxon_file{$1} = $the_filename; # glob($2); # ~/xxx -> /home/tomfy/xxx
#             # print STDERR "BBB: [$1]  [", $taxon_file{$1}, "]\n\n";
#          } elsif (/\S/) {       # not all whitespace
#             warn "B unexpected line in control_file: $_";
#          }
#       } elsif (/^\s*group\s+(\S+)\s+\'([^']*)\'/) {
#          # 'groupname1:a,b,c,d,e;groupname2:h,i,j,k,l'
#          #  print STDERR "QWERTY: [$1] [$2] \n";
#          $added_groups_string .= $1 . ':' . $2 . ';';
#       } elsif (/^\s*(\S+)\s+(\S+)/) { # if param value is not present in file (blank) nothing is stored
#          #  print "param and val: $1 $2 \n";
#          $param_name_val{$1} = $2;
#       } elsif (/\S/) {
#          warn "C unexpected line in control_file: $_";
#       }
#    }
#    #   print STDERR 'Added groups String: ' . "$added_groups_string \n";
#    $added_groups_string =~ s/;\s*$//; # remove final ;
#    $added_groups_string = "'" . $added_groups_string . "'";
#    $param_name_val{added_groups_string} = $added_groups_string;
#    #print STDERR 'Added groups String: ' . "$added_groups_string \n";
#    return (\%taxon_file, \%param_name_val);
# }

sub make_gg{ # make gene-genome string specifying association between seq. ids, and their species.
   my $taxon_file = shift; # hashref; e.g. {'Arabidopsis_lyrata' => 'Alyrata.fasta'}
   my $gg_filename = shift || undef;
   my $gg_string = '';
   my %species_seqcount = ();
   while ( my ($species, $file) = each %$taxon_file) {
      $gg_string .= "$species:";
      my @expanded_filenames = glob $file;
      warn "Using only first expansion: $expanded_filenames[0], of $file.\n" if(scalar @expanded_filenames != 1);
      my $xfile = $expanded_filenames[0];

# print "is [$xfile] a file?  ", (-f $xfile)? 'yes' : 'no', "\n";
# my @fstats = stat $xfile;
# printf("permissions: %o \n", $fstats[2]);
      open my $fh, "<", "$xfile" or die "In make_gg, couldnt open file [$xfile] for reading.\n";
      while (<$fh>) {
         if (/^>\s*(\S+)/) {	# this is a line with a sequence id.
            $gg_string .= " $1";
            $species_seqcount{$species}++;
         }
      }
      close $fh;
      $gg_string .= "\n";
   }
   
   if (defined $gg_filename) {  # write gg file.
      if (open my $fh, ">", "$gg_filename") { 
         print $fh $gg_string; close $fh;
      } else {
         warn "In make_gg. Couldnt open file $gg_filename for writing.\n";
      }
   }
   return \%species_seqcount;
}

sub split_fasta{
   my $fasta_filename = shift;
   my $n_parts = shift;

   #   my ($v, $dir, $fname) = File::Spec->splitpath($fasta_filename);
   #   print STDERR "[$v], [$dir], [$fname] \n";
   my $wcout = `grep '^>' $fasta_filename | wc`;
   $wcout =~ /^\s*(\S+)/;
   my $n_seqs = $1;
   my $n_sequences_in_each_part = int($n_seqs/$n_parts) + 1;

   open my $fh, "<", "$fasta_filename" or warn "couldn't open $fasta_filename for reading.\n";

   my $target_sequence_count = $n_sequences_in_each_part;
   print "total number of sequences: $n_seqs ; number of parts: $n_parts approx. number of sequences per part: $target_sequence_count \n";

   my $i_part = 1;
   my $output_filename_stem = $fasta_filename;
   $output_filename_stem =~ s/[.]fasta$//;
   my $output_filename = $output_filename_stem . "_part$i_part.fasta";
   my @output_fasta_part_filenames = ($output_filename);
   open my $fh_out, ">", $output_filename;
   my $sequence_count = 0;
   while (<$fh>) {
      if (/^>/) {
         s/^(>\s*\S+).*/$1/;    # just keep > and id.
         $sequence_count++;
         #       print "$sequence_count  $_\n";
         if ($sequence_count > $target_sequence_count) {
            close $fh_out;
            $i_part++;
            $output_filename = $output_filename_stem . "_part$i_part.fasta";
            $target_sequence_count += $n_sequences_in_each_part;
            push @output_fasta_part_filenames, $output_filename;
            open $fh_out, ">", $output_filename;
         }
      }
      print $fh_out  $_;
   }
   close $fh_out;
   close $fh;
   return \@output_fasta_part_filenames;
}


sub gg_summary_string{ # returns a string summarizing the gene-genome association information
   my $species_seqcount = shift;
   my $string = "Gene-genome association file $gg_filename created for " . scalar keys (%$species_seqcount) . " species:\n";
   $string .= sprintf("%30s  %8s\n", "species", "   count");
   my ($count_species, $count_seqs) = (0, 0);
   while (my ($sp, $count) = each %$species_seqcount) {
      $string .= sprintf("%30s  %8i\n", $sp, $count);
      $count_species++; $count_seqs += $count;
   }
   $string .= sprintf("totals: %22i  %8i\n", $count_species, $count_seqs);
   return $string;
}

sub get_1sequence_fasta_string{
   my $query_id = shift;
   my $fasta_file = shift;      # , $pl_obj->{query_inputpath});
   my $fasta_string = '';
   my $append_this = 0;
   open my $FHfin, "<", "$fasta_file";
   while (<$FHfin>) {
      if (/^>(\S+)\s/) {
         if ($1 eq $query_id) {
            #     $fasta_string = $_;
            $append_this = 1;
         } else {
            last if($append_this == 1);
            $append_this = 0;
         }
      }
      $fasta_string .= $_ if($append_this);
   }
   return $fasta_string;
}

sub do_this_section{
   my $this_section_name = shift;
   my $start_section_name = shift;
   my $last_section_name = shift;
   my $sectionname_number = shift;
   if ($sectionname_number->{$this_section_name} < $sectionname_number->{$start_section_name}) {
      return 0;
   } elsif ($sectionname_number->{$this_section_name} > $sectionname_number->{$last_section_name}) {
      return 0;
   } else {
      return 1;
   }
}
