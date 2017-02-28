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

my @default_parameters = ( # list of key / value pairs, to initialize a hash  
                          # probably just leave this as empty; defaults are in constructors of objects.
                          # blast_max_matches => 2500, blast_max_e_value => 1e-6,
                    #        min_sequence_length => 12, n_pieces => 2,
                          #  blast_output_dir => 'blast_out',
                          # families_output_dir => 'families',
                          # family_size_limit => '500',
                          # family_max_e_value => 1e-6,
                          # family_multiplicity_knee => 1000,
                          # family_log10_eval_penalty => 0,
                          # family_fastas_output_dir => 'family_fastas',
                          # alignment_output_dir => 'alignments',
                          # alignment_program => 'mafft',
                          # alignment_quality => 'quick',
                          # trees_output_dir => 'newick_trees',
                          # min_nongap_fraction => 0.7,
                          # n_multi_ft => 1,
                          # segments => {'init' => 1, 'blast' => 1, 'make_families' => 1, 'align' => 1},
                         );

my $control_filename = undef;

# Process long cl options
GetOptions(
	   'control_filename=s' => \$control_filename, # required on command line
           #           'start=s' => \$start, # 'blast', 'make_family', etc. 
           #           'last=s' => \$last,
           #           'n_pieces=i' => \$n_pieces,
           #	   'max_eval=s' => \$max_eval,
           #	   'fam_size_limit=i' => \$max_fam_size,
           #           'n_multi_ft=i' => \$n_multi_ft,
           #           'query_id=s' => \$query_id,
           #           'query_taxon=s' => \$cl_query_taxon, # if specified on command line, overrides control file
	  );

die "Control filename must be specified on command line.\n" if(!defined $control_filename);
die "Specified control filename is: $control_filename; file does not exist.\n" if(! -f $control_filename);
my $pl_obj = PhylogenomicPipeline->new($control_filename, \@default_parameters);
print "\n\n\n\n";
print $pl_obj->stringify(), "\n\n";


my $segments = $pl_obj->get('segments');
print join("  ", $segments->keys()), "\n"; # exit;
my $iter = $segments->iterator();
while (my ($segment_name, $segment_obj) = $iter->()) {
   print STDERR  "\n\n\n", "segment name: $segment_name  \n", $segment_obj->stringify(), "\n\n\n\n\n";
   print STDERR "about to call run method of ", $segment_name, " ", $segment_obj, "\n";
   $segment_obj->run();
   print STDERR "after run method of ", $segment_name, " ", $segment_obj, " wd is: ", abs_path('./'),  "\n";
  # print STDERR $segment_obj->stringify();
   #exit;
}
exit;
# ##########################################################################################
# if(1){ # (do_this_section('align_and_make_trees', $start, $last, \%section_number)) {
#    # fork processes to do alignment, tree finding.
#    my @alignment_programs = ('muscle', 'mafft'); # i.e. default is do both
#    if ($alignment_program eq 'muscle') {
#       @alignment_programs = ('muscle');
#    } elsif ($alignment_program eq 'mafft') {
#       @alignment_programs = ('mafft');
#    }

#    for my $align_program (@alignment_programs) {
#       my $align_dir = $align_program . "_" . $pl_obj->{alignments_dir};
#       mkdir $align_dir unless(-d $align_dir);
#       die "Couldnt make directory $align_dir.\n" unless(-d $align_dir); # make a directory to put alignments in.
#       chdir $align_dir or die "Couldnt change directory to $align_dir.\n";
#       my $newick_dir = $pl_obj->{newick_trees_dir};
#       mkdir $newick_dir unless(-d $newick_dir);
#       my $n_alignment_files_to_do = scalar @fastas_filenames;
#       my @alfastas_filenames = ();
#       my @newicks_filenames = ();
#       for my $a_fastas_filename (@fastas_filenames) {
#          my $malign_out_filename = $a_fastas_filename;
#          my ($v, $dir, $fname) = File::Spec->splitpath($malign_out_filename);
#          $malign_out_filename = $fname;
#          $malign_out_filename =~ s/fastas$/alfastas/;
#          push @alfastas_filenames, $malign_out_filename;
#          my $output_newick_filename = $newick_dir . "/" . $malign_out_filename;
#          $output_newick_filename =~ s/alfastas$/newicks/;
#          push @newicks_filenames, $output_newick_filename;
#          my $pid = fork(); # returns 0 to child process, pid of child to parent process.
#          if ($pid == 0) {  # child process
#             $a_fastas_filename = '../' . $a_fastas_filename;
#             $gg_filename = '../' . $gg_filename;
#             print "just before malign. cwd: ", getcwd, "\n";
#             my $malign_cl = "malign.pl  -input $a_fastas_filename  -align $align_program  -quality $alignment_quality  -output $malign_out_filename ";
#             print STDERR "malign command line: [$malign_cl] \n";
#             my $malign_stdout = `$malign_cl`;
#             print "malign.pl finished aligning $a_fastas_filename; output file: $malign_out_filename. \n";

#             my $tree_construct_cl =  "tree_construct.pl -gg $gg_filename -input $malign_out_filename -output $output_newick_filename -nongap_frac $min_nongap_fraction";
#             # "NJFTPHtree.pl -gg $gg_filename -input $malign_out_filename -output $output_newick_filename -n_bs $n_bs -nongap_frac $min_nongap_fraction ";
#             print "tree construction cl: $tree_construct_cl \n";
#             my $tree_construct_stdout = `$tree_construct_cl`;
#             exit(0);
#          }
#       }
#       my $children_done = 0;
#       while (wait() != -1) {
#          $children_done++; print "Number of files finished aligning with $align_program: $children_done out of $n_alignment_files_to_do. \n";
#       }

#       my @files = split(" ", `ls *tmpfile`); # cleanup - delete temp files.
#       for (@files) {
#          unlink $_;
#       }

#       print "cwd: ", getcwd, "\n";
#       chdir '../' or die "chdir ../ failed \n";
#       print "cwd: ", getcwd, "\n";

#       print $fh_progress "Done with aligning with $align_program and constructing trees. \n";
#       print $fh_progress "alfastas files: \n" , join("\n", @alfastas_filenames), "\n";
#       print $fh_progress "newick files: \n", join("\n", @newicks_filenames), "\n";
#       print STDERR "Done with aligning with $align_program and constructing trees. \n";
#       print STDERR "alfastas files: \n" , join("\n", @alfastas_filenames), "\n";
#       print STDERR "newick files: \n", join("\n", @newicks_filenames), "\n";
#    }                            # end of loop over alignment programs
# }                               # end of align_and_make_trees section
# close $fh_progress;
# ##########################################################################################
