#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename 'dirname';
use File::Spec qw(splitpath);
use Cwd qw(abs_path getcwd);
use Time::Piece;

my $control_filename = undef;
my $max_eval = 1e-6;		# default.
my $max_fam_size = 140;
# my $gg_filename_cl = undef;
# my $m8_filename = undef;
my $fasta_infile = undef;
my $n_pieces = 2;
# my $abc_filename = undef;

# my $best_model_only = 0;
my $blast_out_dir = 'blast_out';
my $fams_abc_dir = 'fams_abc';
my $fams_fastas_dir = 'fams_fastas';
my $mafft_alignment_dir = 'mafft_aligment_fastas';
my $muscle_alignment_dir = 'muscle_alignment_fastas';
my $n_multi_ft = 1;
my $cl_query_taxon = undef;
my $query_id = undef;
my $default_min_sequence_length = 20;
my $query_number;               # either 'single' or 'multiple'
# my $added_groups_string = '';

# define some abbreviations for species names, for use in file names
my %species_long_short = ('Medicago_truncatula' => 'Med.tr',
			  'Arabidopsis_thaliana' => 'A.th',
			  'Solanum_lycopersicum' => 'Sol.lyc',
			  'Solanum_tuberosum' => 'Sol.tub',
			  'Brachypodium_distachyon' => 'Brachy.dist',
			  'Sorghum_bicolor' => 'Sorg.bi',
			  'Oryza_sativa' => 'Oryza.sat',
			  'Selaginella_moellendorffii' => 'Sel.mo',
			 );

# Process long cl options
GetOptions(
	   'control_filename=s' => \$control_filename,
	   'max_eval=s' => \$max_eval,
	   'fam_size_limit=i' => \$max_fam_size,
	   'fasta_infile=s' => \$fasta_infile,
	   'n_pieces=i' => \$n_pieces,
           'n_multi_ft=i' => \$n_multi_ft,
           'query_id=s' => \$query_id,
           'query_taxon=s' => \$cl_query_taxon, # if specified on command line, overrides control file
           
	  );


my ($taxon_file, $param_name_val) = get_params_from_control_file($control_filename);
if (defined $cl_query_taxon) {
   $param_name_val->{query_taxon} = $cl_query_taxon;
   $param_name_val->{query_taxon_file} = $taxon_file->{$param_name_val->{query_taxon}};
}

print STDERR "species, file: \n";
while (my ($s, $f) = each %$taxon_file) {
   printf STDERR ("%30s     %50s \n", $s, $f);
}
print "\n";
print STDERR "param, value: \n";
while (my ($p, $v) = each %$param_name_val) {
   #  print STDERR "$p, $v \n";
   printf STDERR ("%30s     %50s \n", $p, $v);

}
print "\n";
#exit;

# get date for incorporating into file names...
my $ltobj = localtime;
my $date = $ltobj->strftime('%Y%h%d');
my $time = $ltobj->hms;
#print STDERR "Pipeline start date, time: $date, ", $ltobj->hms, "\n\n";
print STDERR "Pipeline start date, time: $date, ", $ltobj->hms, "\n";
#exit;

my $qspecies = $param_name_val->{query_taxon};
# $qspecies = $species_long_short{$qspecies} if(exists $species_long_short{$qspecies});
print STDERR "qspecies: $qspecies \n";
my $n_species = scalar keys %$taxon_file;
my $filename_head = "$qspecies" . "_vs_" . "$n_species" . "_" . $date;
print STDERR "filename head: $filename_head \n";
my $progress_filename = $filename_head . '.progress';
open my $fh_progress, ">", "$progress_filename";



my $blast_max_matches = (exists $param_name_val->{blast_max_matches})? $param_name_val->{blast_max_matches} : 2500;
my $blast_max_e_value = (exists $param_name_val->{blast_max_e_value})? $param_name_val->{blast_max_e_value} : '1e-8';

my $family_size_limit = (exists $param_name_val->{family_size_limit})? $param_name_val->{family_size_limit} : undef;
my $family_max_e_value = (exists $param_name_val->{family_max_e_value})? $param_name_val->{family_max_e_value} : '1e-12';
my $family_multiplicity_knee =  (exists $param_name_val->{family_multiplicity_knee})? $param_name_val->{family_multiplicity_knee} : 6;
my $family_log10_eval_penalty = (exists $param_name_val->{family_log10_eval_penalty})? $param_name_val->{family_log10_eval_penalty} : 12;
my $family_taxon_requirement = (exists $param_name_val->{family_taxon_requirement})? 
  $param_name_val->{family_taxon_requirement} : '7dicots,6; 4monocots,3';
my $added_groups_string = (exists $param_name_val->{added_groups_string})?  $param_name_val->{added_groups_string} : undef;

my $alignment_program = (exists $param_name_val->{alignment_program})? $param_name_val->{alignment_program} : 'both'; # muscle, mafft, or both.
my $alignment_quality = (exists $param_name_val->{alignment_quality})? $param_name_val->{alignment_quality} : 'best'; # quick or best.

my $min_nongap_fraction = (exists $param_name_val->{min_nongap_fraction})? $param_name_val->{min_nongap_fraction} : 0.15;

$n_multi_ft = $param_name_val->{n_multi_ft} if(defined $param_name_val->{n_multi_ft});
$n_multi_ft = ($n_multi_ft > 1)? $n_multi_ft : 1;
my $n_bs = ($n_multi_ft >= 1)? $n_multi_ft - 1 : 0;
# ********** make gg string and file (gene-genome association)
my $gg_filename = "$filename_head.gg";
my ($gg_string, $species_seqcount) = make_gg($taxon_file, $gg_filename);
my $gg_summary = gg_summary_string($species_seqcount);
print STDERR $gg_summary;
print $fh_progress $gg_summary;
print $fh_progress "gg file: $gg_filename \n";

# ********** make blast db ****************
my $all_species_fasta_filename = "all_" . $n_species . "_species_" . $date . ".fasta";
my $fasta_files = '';
for my $filename (values %$taxon_file) {
   $fasta_files .= "$filename ";
}

my $asff = $all_species_fasta_filename . "_x";
system "cat $fasta_files | clean_fasta_idlines.pl > $asff";
my $min_seq_length = (defined $param_name_val->{min_sequence_length})? $param_name_val->{min_sequence_length} :  $default_min_sequence_length;
print STDERR "About to call remove_short_seqs.pl with min_seq_length = $min_seq_length \n";
system "remove_short_seqs.pl $min_seq_length < $asff > $all_species_fasta_filename";
system "formatdb -p T -i $all_species_fasta_filename ";
print STDERR "blast db created for $all_species_fasta_filename.\n";
print $fh_progress "blast db created for $all_species_fasta_filename.\n";

my $qfasta_filename;
my $fasta_part_filenames;
# if doing single sequence:
if (defined $query_id) {
   $query_number = 'single';
   my $query_fasta_string = get_1sequence_fasta_string($query_id, $param_name_val->{query_inputpath});
   print STDERR "AAAAAAAAAAAAa: [[$query_fasta_string]]\n";
   my @lines = split("\n", $query_fasta_string);
   my $id_line = shift @lines;
   my $id = ($id_line =~ /^>(\S+)/)? $1 : 'XXX_';
   my $query_filename = $id . ".fasta";
   $fasta_part_filenames = [$query_filename];
   open my $FHqff, ">", "$query_filename";
   $query_fasta_string = $id_line . "\n" . join("", @lines) . "\n";
   print $FHqff $query_fasta_string, "\n";
} else {
   # else doing whole genome ('all'):
   $query_number = 'multiple';
   # ********* split query fasta **************
   $qfasta_filename = $param_name_val->{query_inputpath};
   my ($vol, $path, $qfname) = File::Spec->splitpath($qfasta_filename);
   #print "[$vol]  [$path]  [$fname]   \n";
 #  my $qfname = $qfasta_filename;
   $qfname =~ s/[.]fasta//;   # remove final .fasta
   #$qfname =~ $param_name_val->{} . '/' . 
   my $qfasta_filename_noshortseqs = $qfname . "_nss.fasta"; # just goes in cwd
   system "remove_short_seqs.pl $min_seq_length < $qfasta_filename > $qfasta_filename_noshortseqs"; # remove short sequences 
   $fasta_part_filenames = split_fasta($qfasta_filename_noshortseqs, $param_name_val->{n_pieces});
   print STDERR "query part fasta files: ", join(", ", @$fasta_part_filenames), "\n";
}


# ********** run blast *********************
my @blast_out_m8_filenames = ();
for my $q_fasta_part_filename (@$fasta_part_filenames) { # loop over the parts which the set of queries has been divided into.
   my $blast_out_filename = $q_fasta_part_filename;
   $blast_out_filename =~ s/fasta$/m8/;
   my ($v, $dir, $fname) = File::Spec->splitpath($blast_out_filename);
   my $blastout_dir = $param_name_val->{blastout_dir};
   mkdir $blastout_dir unless(-d $blastout_dir);
   $blast_out_filename = $blastout_dir . '/' . $fname;
   print "blastout path: $blast_out_filename \n";
   push @blast_out_m8_filenames, $blast_out_filename;
   my $pid = fork(); # returns 0 to child process, pid of child to parent process.
   if ($pid == 0) {  # child process
   
      # blastall -p blastp -d dbfilename -i queryfilename -e ‘1e-6’ -m 8 -a 4 >  xx_blastout.m8
      my $blast_cl = "nohup blastall -p blastp -d $all_species_fasta_filename -i $q_fasta_part_filename -e '$blast_max_e_value' -b $blast_max_matches -m 8 -o $blast_out_filename";
      print "blast cl: $blast_cl \n";
      my $blast_stdout = `$blast_cl`;
      print "blast finished analyzing $q_fasta_part_filename; output file: $blast_out_filename. \n";
      exit(0);
   }
}				# end loop over parts
my $children_done = 0;
while (wait() != -1) {
   $children_done++; print "Number of files finished analyzing with blast: $children_done. \n";
}
print STDERR "blast finished. blast output filenames: \n", join("\n", @blast_out_m8_filenames), "\n";
print $fh_progress "blast finished. blast output filenames: \n", join("\n", @blast_out_m8_filenames), "\n";

#exit;
# *************** done running blast ******************

my @abc_part_filenames = ();
for my $m8_filename (@blast_out_m8_filenames) {
   my ($v, $dir, $fname) = File::Spec->splitpath($m8_filename);
   $fname =~ s/m8$/abc/;
   my $fam_abcs_dir = $param_name_val->{fam_abcs_dir};
   mkdir $fam_abcs_dir unless(-d $fam_abcs_dir); # create dir if doesn't exist.
   my $fam_abcs_filename = $fam_abcs_dir . '/' . $fname;
   my $m8tofams_cl = "m8toabc_fams.pl -input $m8_filename -max_eval $family_max_e_value -gg $gg_filename -output $fam_abcs_filename ";
   $m8tofams_cl .= "-multiplicity_knee $family_multiplicity_knee -log10penalty $family_log10_eval_penalty ";
   my $abc_filename = `$m8tofams_cl`;
   # -fam_size_limit $family_size_limit `;
   $abc_filename =~ s/\s+$//; 
   print STDERR "[$abc_filename]\n";
   push @abc_part_filenames, $abc_filename;
}
#my $abc_part_filenames = `split_abc.pl $abc_filename $n_pieces`;
#my @abc_parts = split(" ", $abc_part_filenames);
#exit;

print "[" ,join("][", @abc_part_filenames), "]\n";
print $fh_progress "abc family filenames: \n", join("\n", @abc_part_filenames), "\n";

# ************** done making abc files *******************
# exit;
my @fastas_filenames = ();
die "fasta input file undefined or file does not exist.\n" unless(defined $all_species_fasta_filename and -f $all_species_fasta_filename);
my $fam_fastas_dir = $param_name_val->{fam_fastas_dir};
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
# exit;

# fork processes to do alignment, tree finding.
my @alignment_programs = ('muscle', 'mafft'); # i.e. default is do both
if ($alignment_program eq 'muscle') {
   @alignment_programs = ('muscle');
} elsif ($alignment_program eq 'mafft') {
   @alignment_programs = ('mafft');
}

for my $align_program (@alignment_programs) {
   my $align_dir = $align_program . "_" . $param_name_val->{alignments_dir};
   mkdir $align_dir unless(-d $align_dir);
   die "Couldnt make directory $align_dir.\n" unless(-d $align_dir); # make a directory to put alignments in.
   chdir $align_dir or die "Couldnt change directory to $align_dir.\n";
   my $newick_dir = $param_name_val->{newick_trees_dir};
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
}                               # end of loop over alignment programs
close $fh_progress;


sub get_params_from_control_file{
   my $control_filename = shift;
   open my $fh_ctrl, "<", "$control_filename" or die "Could not open $control_filename for reading; exiting.";
   my $added_groups_string = '';
   my %param_name_val = ();
   my %taxon_file = ();
   while (<$fh_ctrl>) {
      next if(/^\s*#/);         # skip all comment lines
      if (/^\s*query_taxon_inputpath/) {
         if (/^\s*query_taxon_inputpath\s+(\S+)\s+(\S+)/) { # $1 is taxon name, e.g. Zea_mays, $2 is path to corresponding input file (protein-models fasta format)
            #    $taxon_file{$1} = $2;
            $param_name_val{query_taxon} = $1;
            $param_name_val{query_inputpath} = $2;
         } elsif (/\S/) {       # not all whitespace
            warn "A unexpected line in control_file: $_";
         }
      } elsif (/^\s*taxon_inputpath/) {
         if (/^\s*taxon_inputpath\s+(\S+)\s+(\S+)/) {
            my @exp_filenames = glob($2);
            my $the_filename = $exp_filenames[0]; 
            #   print STDERR "AAA: [", $2, "]  [", glob($2), "] [$the_filename]\n";            
            $taxon_file{$1} = $the_filename; # glob($2); # ~/xxx -> /home/tomfy/xxx
            # print STDERR "BBB: [$1]  [", $taxon_file{$1}, "]\n\n";
         } elsif (/\S/) {       # not all whitespace
            warn "B unexpected line in control_file: $_";
         }
      } elsif (/^\s*group\s+(\S+)\s+\'([^']*)\'/) {
         # 'groupname1:a,b,c,d,e;groupname2:h,i,j,k,l'
       #  print STDERR "QWERTY: [$1] [$2] \n";
         $added_groups_string .= $1 . ':' . $2 . ';';
      } elsif (/^\s*(\S+)\s+(\S+)/) { # if param value is not present in file (blank) nothing is stored
         #  print "param and val: $1 $2 \n";
         $param_name_val{$1} = $2;
      } elsif (/\S/) {
         warn "C unexpected line in control_file: $_";
      }
   }
#   print STDERR 'Added groups String: ' . "$added_groups_string \n";
   $added_groups_string =~ s/;\s*$//; # remove final ;
$added_groups_string = "'" . $added_groups_string . "'";
   $param_name_val{added_groups_string} = $added_groups_string;
#print STDERR 'Added groups String: ' . "$added_groups_string \n";
   return (\%taxon_file, \%param_name_val);
}

sub make_gg{ # make gene-genome string specifying association between seq. ids, and their species.
   my $taxon_file = shift;	# hashref 
   my $gg_filename = shift || undef;
   my $gg_string = '';
   my %species_seqcount = ();
   while ( my ($species, $file) = each %$taxon_file) {
      $gg_string .= "$species:";
      open my $fh, "<", "$file" or die "In make_gg, couldnt open file [$file] for reading.\n";
      while (<$fh>) {
         if (/^>\s*(\S+)/) {	# this is a line with a sequence id.
            $gg_string .= " $1";
            $species_seqcount{$species}++;
         }
      }
      $gg_string .= "\n";
   }
   if (defined $gg_filename) {  # write gg file.
      if (open my $fh, ">", "$gg_filename") { 
         print $fh $gg_string; close $fh;
      } else {
         warn "In make_gg. Couldnt open file $gg_filename for writing.\n";
      }
   }
   return ($gg_string, \%species_seqcount);
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
   my $fasta_file = shift;    # , $param_name_val->{query_inputpath});
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
