package PipeSegBlast;
use strict;
use base 'PipelineSegment';
use List::Util qw (min max sum);
use TomfyMisc qw (date_time split_fasta);

sub new {
   my $class = shift;
   my $pl_obj = shift;          # 
   my $args  = {};
   my $self  = bless $args, $class;
   $self->set('pipeline', $pl_obj); # store the pipeline obj to which this segment belongs.
   return $self;
}

sub run{
   my $self = shift;
   # **************** concatenate all target fasta files, 'clean' id lines, and remove short sequences *************  
   my $fasta_files = '';
   for my $filename (values %{$self->get('pipeline')->get('taxon_inputpath')}) {
      $fasta_files .= "$filename ";
   }

   my ($date, $time) = date_time();
my $min_seq_length = $self->get_plattr('min_sequence_length');
my $fh_progress = $self->get_plattr('fh_prog');
   my $all_species_fasta_filename = "all_" . $self->get('pipeline')->get('n_species') . "_species_" . $date . ".fasta";
   my $asff = $all_species_fasta_filename . "_x";
   system "cat $fasta_files | clean_fasta_idlines.pl > $asff";
   print STDERR "About to call remove_short_seqs.pl with min_seq_length = $min_seq_length \n";
   system "remove_short_seqs.pl $min_seq_length < $asff > $all_species_fasta_filename";

   # ********** make blast db ****************
   system "formatdb -p T -i $all_species_fasta_filename ";
   print STDERR "blast db created for $all_species_fasta_filename.\n";
   print $fh_progress "blast db created for $all_species_fasta_filename.\n";
   # ******** Done making blast db ***************


#   my $qfasta_filename;
   my $fasta_part_filenames;
   # if doing single sequence:
   my $query_id = undef;
   if (defined $query_id) {
      die "single query option not implemented.\n";
      # $query_number = 'single';
      # my $query_fasta_string = get_1sequence_fasta_string($query_id, $pl_obj->get('query_inputpath'));
      # print STDERR "AAAAAAAAAAAAa: [[$query_fasta_string]]\n";
      # my @lines = split("\n", $query_fasta_string);
      # my $id_line = shift @lines;
      # my $id = ($id_line =~ /^>(\S+)/)? $1 : 'XXX_';
      # my $query_filename = $id . ".fasta";
      # $fasta_part_filenames = [$query_filename];
      # open my $FHqff, ">", "$query_filename";
      # $query_fasta_string = $id_line . "\n" . join("", @lines) . "\n";
      # print $FHqff $query_fasta_string, "\n";
   } else {
      # else doing whole genome ('all'):
    #  $query_number = 'multiple';
      # ********* split query fasta **************
      my @qts = ();
      my @qtpaths = ();
      while (my($qt, $qtpath) = each %{$self->get_plattr('query_taxon_inputpath')}) {
         push @qts, $qt;
         push @qtpaths, $qtpath;
      }
      my $qfasta_filename = join(".", @qts) . ".fasta";
      my $qtpaths_str = join(" ", @qtpaths);
      system "cat $qtpaths_str > $qfasta_filename";
      $self->set('query_inputpath'); # not needed?
      my ($vol, $path, $qfname) = File::Spec->splitpath($qfasta_filename);
      $qfname =~ s/[.]fasta//;  # remove final .fasta
      my $qfasta_filename_noshortseqs = $qfname . "_nss.fasta"; # just goes in cwd
      system "remove_short_seqs.pl $min_seq_length < $qfasta_filename > $qfasta_filename_noshortseqs"; # remove short sequences 
      $fasta_part_filenames = split_fasta($qfasta_filename_noshortseqs, $self->get_plattr('n_pieces'));
      print STDERR "query part fasta files: ", join(", ", @$fasta_part_filenames), "\n";
   }

my @blast_out_m8_filenames = ();
   # ********** run blast *********************
   for my $q_fasta_part_filename (@$fasta_part_filenames) { # loop over the parts which the set of queries has been divided into.
      my $blast_out_filename = $q_fasta_part_filename;
      $blast_out_filename =~ s/fasta$/m8/;
      my ($v, $dir, $fname) = File::Spec->splitpath($blast_out_filename);
      my $blastout_dir = $self->get_plattr('blast_output_dir');
      mkdir $blastout_dir unless(-d $blastout_dir);
      $blast_out_filename = $blastout_dir . '/' . $fname;
      print "blastout path: $blast_out_filename \n";
      push @blast_out_m8_filenames, $blast_out_filename;
      my $pid = fork(); # returns 0 to child process, pid of child to parent process.
      if ($pid == 0) {  # child process
   
         # blastall -p blastp -d dbfilename -i queryfilename -e ‘1e-6’ -m 8 -a 4 >  xx_blastout.m8
         my $blast_max_matches = $self->get_plattr('blast_max_matches');
         my $blast_max_e_value = $self->get_plattr('blast_max_e_value');
         my $blast_cl = "nohup blastall -p blastp -d $all_species_fasta_filename -i $q_fasta_part_filename -e $blast_max_e_value -b $blast_max_matches -m 8 -o $blast_out_filename";
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

}



1;
