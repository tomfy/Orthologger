package PipeSegBlast;
use strict;
use base 'PipelineSegment';
use List::Util qw (min max sum);
use TomfyMisc qw (date_time split_fasta_file split_fasta_string fix_fasta_files);
use Cwd qw (getcwd abs_path);


sub new {
   my $class = shift;
   my $pl_obj = shift;
   print STDERR "top of pipesegblast constructor.\n";
   my $self  = $class->SUPER::new($pl_obj); # bless $args, $class;
   $self->init({
                output_dir => 'blast_out',
                max_matches => 2500,
                max_e_value => 1e-6,
               });
   print STDERR "bottom of pipesegblast constructor.\n";
   return $self;
}



sub run{
   my $self = shift;

   my $pl_obj = $self->get('pipeline');
   my ($date, $time) = date_time();
   my $min_seq_length = $self->get('min_sequence_length');
   my $fh_progress = $pl_obj->get('fh_progress');
   my $all_species_fasta_filename = abs_path( "all_" . $pl_obj->get('n_species') . "_species_" . $date . ".fasta" );

   # **************** concatenate all target fasta files, 'clean' id lines, and remove short sequences *************  
   my @fasta_files = $pl_obj->get('taxon_inputpath')->values();
   open my $fhout, ">", "$all_species_fasta_filename" or die "Couldn't open $all_species_fasta_filename for writing.\n";
   print $fhout fix_fasta_files($min_seq_length, @fasta_files);
   close $fhout;

   # ********** make blast db ****************
   system "formatdb -p T -i $all_species_fasta_filename ";
   print STDERR "blast db created for $all_species_fasta_filename.\n";
   print $fh_progress "blast db created for $all_species_fasta_filename.\n";
   # ******** Done making blast db ***************

   # ********* split query fasta **************
   my @q_fasta_part_filenames = ();
   my $query_fasta_string = fix_fasta_files($min_seq_length, $pl_obj->get('querytaxon_inputpath')->values());
   my $q_fasta_strings = split_fasta_string($query_fasta_string, $pl_obj->get('n_pieces'));
   while (my ($i, $v) = each @$q_fasta_strings) {
      my $part_filename = abs_path(join('', @{$pl_obj->get('querytaxa_short')}) . "_" . $date . "_" . 'part' . ($i+1) . ".fasta");
      push @q_fasta_part_filenames, $part_filename;
      open $fhout, ">", "$part_filename" or die "Couldn't open $part_filename for writing.\n";
      print $fhout $v;
      close $fhout;
   }
#print join("\n", @q_fasta_part_filenames), "\n"; exit;

   my $blastout_dir = $self->get('output_dir');
   mkdir $blastout_dir unless(-d $blastout_dir);
chdir($blastout_dir) or die "Couldn't chdir to $blastout_dir.\n";
   my @blast_out_m8_filenames = ();
   # ********** run blast *********************
   for my $q_fasta_part_filename (@q_fasta_part_filenames) { # loop over the parts which the set of queries has been divided into.
      my ($v, $dir, $blast_out_filename) = File::Spec->splitpath($q_fasta_part_filename);
      $blast_out_filename =~ s/^[^_]+//;
      $blast_out_filename =~ s/[.]fasta/.m8/;
      $blast_out_filename = 
      join('', @{$pl_obj->get('querytaxa_short')}) . "_v_" . $pl_obj->get('n_pieces') . $blast_out_filename; # "_" . $date . 'part' . "_" . ($i+1) . ".fasta"
      $blast_out_filename = abs_path($blast_out_filename);
 print "$q_fasta_part_filename $all_species_fasta_filename \n";
      print "blastout path: $blast_out_filename \n";
      push @blast_out_m8_filenames, $blast_out_filename;
      my $pid = fork(); # returns 0 to child process, pid of child to parent process.
      if ($pid == 0) {  # child process
   
         # blastall -p blastp -d dbfilename -i queryfilename -e ‘1e-6’ -m 8 -a 4 >  xx_blastout.m8
         my $blast_max_matches = $self->get('max_matches');
         my $blast_max_e_value = $self->get('max_e_value');
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
   my $blastout_filenames = join(" ", @blast_out_m8_filenames);
   $self->set('blast_out_m8_filenames', $blastout_filenames);
   print STDERR "blast finished. blast output filenames: \n", $self->get('blast_out_m8_filenames'), "\n";
   print $fh_progress "blast finished. blast output filenames: \n", $self->get('blast_out_m8_filenames'), "\n";

 #  chdir($blastout_dir) or die "Couldn't chdir to $blastout_dir.\n";

   my $state_filename = 'pipeline_state';
   open my $fh, ">", "$state_filename" or die "couldn't open $state_filename for writing.\n";
   print $fh $pl_obj->stringify();
   # *************** done running blast ******************

}


1;
