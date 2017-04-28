#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename 'dirname';
use File::Spec qw(splitpath);
use Cwd qw(abs_path getcwd);
use Time::Piece;
use constant STATE_FILENAME => 'pipeline_state';

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

my @default_parameters = (
                          # list of key / value pairs, to initialize a hash  
                          # probably just leave this as empty; defaults are in constructors of objects.
                          state_filename => STATE_FILENAME,
                          # blast_max_matches => 2500, blast_max_e_value => 1e-6,
                          #   min_sequence_length => 12, n_pieces => 2,
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

my $pl_obj;
my $upstream_state_filename = "../" . STATE_FILENAME; # pipeline_state";
if (-f $upstream_state_filename) { # this is the case of restarting from the middle of pipeline
   print "using upstream state file:  '$upstream_state_filename'\n"; #exit;
   $pl_obj = PhylogenomicPipeline->new($upstream_state_filename, @default_parameters);
   $pl_obj->reset_start_output_dir('./');
   # for ($pl_obj->get('segments')->values()) {
   #    print "seg name, state: ", $_->get('segment_name'), "  ", $_->get('state'), "\n";
   # }
} else {
   print "using control file from command line. '$control_filename'\n"; #exit;
   die "Control filename must be specified on command line.\n" if(!defined $control_filename);
   die "Specified control filename is: $control_filename; file does not exist.\n" if(! -f $control_filename);
   $pl_obj = PhylogenomicPipeline->new($control_filename, @default_parameters);
}
print "\n\n";

$pl_obj->run();

# end
