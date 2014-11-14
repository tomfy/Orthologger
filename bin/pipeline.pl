#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Time::Piece;


my $control_filename = undef;
my $max_eval = 1e-12;		# default.
my $max_fam_size = 140;
# my $gg_filename_cl = undef;
# my $m8_filename = undef;
my $fasta_infile = undef;
my $n_pieces = 8;
# my $abc_filename = undef;
my $first_of_species_factor = 1; # has no effect if no ggfile.
# by setting this to some large number (e.g. 1e20) will keep the best match
# from each species, if it is within this factor of the e-value threshold
# so if it is 1, it has no effect.

# my $best_model_only = 0;

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
	   #	   'input_filename=s' => \$m8_filename,
	   'max_eval=s' => \$max_eval,
	   'fam_size_limit=i' => \$max_fam_size,
	   #	   'ggfilename=s' => \$gg_filename_cl,
	   'fasta_infile=s' => \$fasta_infile,
	   # 'only_best_model=s' => \$best_model_only
	   'n_pieces=i' => \$n_pieces,
	  );


my ($file_taxon, $param_name_val) = get_params_from_control_file($control_filename);

print STDERR "species, file: \n";
while (my ($f, $s) = each %$file_taxon) {
  print STDERR "$s, $f\n";
}
print "\n";
print STDERR "param, value: \n";
while (my ($p, $v) = each %$param_name_val) {
  print STDERR "$p, $v \n";
}
print "\n";


# get date for incorporating into file names...
my $ltobj = localtime;
my $date = $ltobj->strftime('%Y%h%d');
my $time = $ltobj->hms;
print STDERR "Pipeline start date, time: $date, ", $ltobj->hms, "\n\n";
print STDERR "Pipeline start date, time: $date, ", $ltobj->hms, "\n";
# exit;

my $qspecies = $file_taxon->{$param_name_val->{query_fasta_filename}};
$qspecies = $species_long_short{$qspecies} if(exists $species_long_short{$qspecies});
print STDERR "qspecies: $qspecies \n";
my $n_species = scalar keys %$file_taxon;
my $filename_head = "$qspecies" . "_vs_" . "$n_species" . "_" . $date;
print STDERR "filename head: $filename_head \n";
my $progress_filename = $filename_head . '.progress';
open my $fh_progress, ">", "$progress_filename";



my $blast_max_matches = (exists $param_name_val->{blast_max_matches})? $param_name_val->{blast_max_matches} : 2500;
my $blast_max_e_value = (exists $param_name_val->{blast_max_e_value})? $param_name_val->{blast_max_e_value} : '1e-8';
my $family_size_limit = (exists $param_name_val->{family_size_limit})? $param_name_val->{family_size_limit} : 140;
my $family_max_e_value = (exists $param_name_val->{family_max_e_value})? $param_name_val->{family_max_e_value} : '1e-12';
my $family_taxon_requirement = (exists $param_name_val->{family_taxon_requirement})? 
  $param_name_val->{family_taxon_requirement} : '7dicots,6; 4monocots,3';
my $alignment_program = (exists $param_name_val->{alignment_program})? $param_name_val->{alignment_program} : 'both'; # muscle, mafft, or both.
my $alignment_quality = (exists $param_name_val->{alignment_quality})? $param_name_val->{alignment_quality} : 'best'; # quick or best.

  # ********** make gg string and file (gene-genome association)
  my $gg_filename = "$filename_head.gg";
my ($gg_string, $species_seqcount) = make_gg($file_taxon, $gg_filename);
my $gg_summary = gg_summary_string($species_seqcount);
print STDERR $gg_summary;
print $fh_progress $gg_summary;
print $fh_progress "gg file: $gg_filename \n";

# ********** make blast db ****************
my $all_species_fasta_filename = "all_" . $n_species . "_species_" . $date . ".fasta";
my $fasta_files = '';
for my $filename (keys %$file_taxon) {
  $fasta_files .= "$filename ";
}

my $asff = $all_species_fasta_filename . "_x";
system "cat $fasta_files | clean_fasta_idlines.pl > $asff";
system "remove_short_seqs.pl < $asff > $all_species_fasta_filename";
system "formatdb -p T -i $all_species_fasta_filename ";
print STDERR "blast db created for $all_species_fasta_filename.\n";
print $fh_progress "blast db created for $all_species_fasta_filename.\n";

# ********* split query fasta **************
my $qfasta_filename = $param_name_val->{query_fasta_filename};
my $qfasta_filename_noshortseqs = $qfasta_filename . "_nss";
system "remove_short_seqs.pl < $qfasta_filename > $qfasta_filename_noshortseqs"; # remove short sequences 
my $fasta_part_filenames = split_fasta($qfasta_filename_noshortseqs, $param_name_val->{n_parts});

# ********** run blast *********************
my @blast_out_m8_filenames = ();
for my $q_fasta_part_filename (@$fasta_part_filenames) { # loop over the parts which the set of queries has been divided into.
  my $blast_out_filename = $q_fasta_part_filename;
  $blast_out_filename =~ s/fasta$/m8/;
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


# *************** done running blast ******************

my @abc_part_filenames = ();
for my $m8_filename (@blast_out_m8_filenames) {
  my $abc_filename = `m8toabc_fams.pl -input $m8_filename -max_eval $family_max_e_value -fam_size_limit $family_size_limit`;
  $abc_filename =~ s/\s+$//; 
  print STDERR "[$abc_filename]\n";
  push @abc_part_filenames, $abc_filename;
}
#my $abc_part_filenames = `split_abc.pl $abc_filename $n_pieces`;
#my @abc_parts = split(" ", $abc_part_filenames);

print "[" ,join("][", @abc_part_filenames), "]\n";
print $fh_progress "abc family filenames: \n", join("\n", @abc_part_filenames), "\n";

# ************** done making abc files *******************
exit;
my @fastas_filenames = ();
die "fasta input file undefined or file does not exist.\n" unless(defined $all_species_fasta_filename and -f $all_species_fasta_filename); 
for my $the_abc_part (@abc_part_filenames) {
  my $output_fastas_filename = $the_abc_part;
  $output_fastas_filename =~ s/abc$/fastas/;
my $seqm2f_cl = "seq+matches2fasta.pl -gg $gg_filename -abc_file $the_abc_part -fasta_infile $all_species_fasta_filename -output_filename $output_fastas_filename";
$seqm2f_cl .= " -taxon_requirement '$family_taxon_requirement' " if(defined $family_taxon_requirement);
  print STDERR "seq+... cl:  $seqm2f_cl \n";
  system("$seqm2f_cl");
# "seq+matches2fasta.pl -gg $gg_filename -abc_file $the_abc_part -fasta_infile $all_species_fasta_filename -output_filename $output_fastas_filename");
  push @fastas_filenames, $output_fastas_filename
}

print STDERR "Family fastas file ready to align: \n", join("\n", @fastas_filenames), "\n";

print $fh_progress "Family_fastas_files: " , join("  ", @fastas_filenames), "\n";
# exit;

# fork processes to do alignment, tree finding.
my @alignment_programs = ('muscle', 'mafft');
if($alignment_program eq 'muscle'){
@alignment_programs = ('muscle');
}elsif($alignment_program eq 'mafft'){
@alignment_programs = ('mafft');
}

for my $align_program (@alignment_programs) {
  mkdir $align_program or die "Couldnt make directory $align_program.\n"; # make a directory to put alignments, etc. in.
  chdir $align_program or die "Couldnt change directory to $align_program.\n";
  
  my $n_alignment_files_to_do = scalar @fastas_filenames;
  my @alfastas_filenames = ();
my @newicks_filenames = ();
  for my $a_fastas_filename (@fastas_filenames) {
    my $malign_out_filename = $a_fastas_filename;

    $malign_out_filename =~ s/fastas$/alfastas/;
    push @alfastas_filenames, $malign_out_filename;
    my $output_newick_filename = $malign_out_filename;
      $output_newick_filename =~ s/alfastas$/newicks/;
    push @newicks_filenames, $output_newick_filename;
    my $pid = fork(); # returns 0 to child process, pid of child to parent process.
    if ($pid == 0) {  # child process
      $a_fastas_filename = '../' . $a_fastas_filename;
      $gg_filename = '../' . $gg_filename;
      my $malign_cl = "malign.pl  -input $a_fastas_filename  -align $align_program  -quality $alignment_quality  -output $malign_out_filename ";
      print STDERR "malign command line: [$malign_cl] \n";
      my $malign_stdout = `$malign_cl`;
      print "malign.pl finished aligning $a_fastas_filename; output file: $malign_out_filename. \n";

      my $nj_ft_bs_stdout = `nj_ft_bs.pl -gg $gg_filename -input $malign_out_filename -output $output_newick_filename`;
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
} # end of loop over alignment programs
close $fh_progress;


sub get_params_from_control_file{
  my $control_filename = shift;
  open my $fh_ctrl, "<", "$control_filename" or die "Could not open $control_filename for reading; exiting.";
  my %param_name_val = ();
  my %file_taxon = ();
  while (<$fh_ctrl>) {
    next if(/^\s*#/); # skip all comment lines
    if (/^\s*query_fasta_filename/) {
      if (/^\s*query_fasta_filename\s+(\S+)\s+(\S+)/) {
	$file_taxon{$1} = $2;
	$param_name_val{query_fasta_filename} = $1;
      } elsif (/\S/) {		# not all whitespace
	warn "A unexpected line in control_file: $_";
      }
    } elsif (/^\s*fasta_filename/) {
      if (/^\s*fasta_filename\s+(\S+)\s+(\S+)/) {
	$file_taxon{$1} = $2;
      } elsif (/\S/) {		# not all whitespace
	warn "B unexpected line in control_file: $_";
      }
    } elsif (/^\s*(\S+)\s+(\S+)/) { # if param value is not present in file (blank) nothing i sstored
   #  print "param and val: $1 $2 \n";
      $param_name_val{$1} = $2; 
    } elsif (/\S/) {
      warn "C unexpected line in control_file: $_";
    }
  }
  return (\%file_taxon, \%param_name_val);
}


sub make_gg{ # make gene-genome string specifying association between seq. ids, and their species.
  my $file_species = shift;	# hashref 
  my $gg_filename = shift || undef;
  my $gg_string = '';
my %species_seqcount = ();
  while ( my ($file, $species) = each %$file_species) {
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
  if (defined $gg_filename) { # write gg file.
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

  my $wcout = `grep '^>' $fasta_filename | wc`;
  $wcout =~ /^\s*(\S+)/;
  my $n_seqs = $1;
  my $n_sequences_in_each_part = int($n_seqs/$n_parts) + 1;

  open my $fh, "<", "$fasta_filename";

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
      s/^(>\s*\S+).*/$1/; # just keep > and id.
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
while(my ($sp, $count) = each %$species_seqcount){
  $string .= sprintf("%30s  %8i\n", $sp, $count);
  $count_species++; $count_seqs += $count;
}
$string .= sprintf("totals: %22i  %8i\n", $count_species, $count_seqs);
return $string;
}
