package PhylogenomicPipeline;
use strict;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use File::Basename 'dirname';
#use File::Spec qw(splitpath);
use Cwd qw(abs_path getcwd);use TomfyMisc qw(date_time short_species_name clean read_block);
use PipeSegGG;
use PipeSegBlast;
use PipeSegRbmcl;
#use PipeSegFamilies;

use parent 'Pipeline';

sub new {
   my $class = shift;
   my $control_filename = shift; #
   my @args = @_;  # e.g. (min_sequence_length => 15, n_pieces => 14);
   print STDERR "top of PhylogenomicsPipeline constructor. control_filename: $control_filename, args:[|", join(";", @args), "|]\n";
   my $default_params = Hash::Ordered->new(
                                           hidden => undef, # OrderedHash->new(), # seqid_species => 1, gg_string => 1, fh_progress => 1),
                                           min_sequence_length => 20,
                                           n_pieces => 3,
                                          );
   $default_params->merge(@args);
   print STDERR "In PhylogenomicPipeline constructor, before calling SUPER::new \n";
   my $self = $class->SUPER::new();
   print STDERR "seqid_species is ", $self->_hidden('seqid_species')? 'HIDDEN' : 'not hidden', "\n"; #exit;
   print STDERR "In PhylogenomicPipeline constructor, after calling SUPER::new \n";
   $self->SUPER::initialize($control_filename, $default_params);
   $self->hide('seqid_species', 'gg_string', 'fh_progress');
   print STDERR "In PhylogenomicPipeline constructor, after calling SUPER::initialize \n";
   print STDERR "seqid_species is ", $self->_hidden('seqid_species')? 'HIDDEN' : 'not hidden', "\n";
   $self->initialize();
   print STDERR "In PhylogenomicPipeline constructor, after calling initialize \n";

   print STDERR "seqid_species is ", $self->_hidden('seqid_species')? 'HIDDEN' : 'not hidden', "\n"; #exit;

   return $self;
}

sub initialize{
   my $self = shift;
   print STDERR "Top of PhylogenomicPipeline::initialize() \n";

   $self->set('n_species',  scalar $self->get('taxon_inputpath')->keys()); # scalar keys %{$self->get('taxon_inputpath')});

   my $qsp_files = $self->get('querytaxon_inputpath') // die "querytaxon_inputpath Hash::Ordered undefined.\n";
   if (blessed $qsp_files and $qsp_files->isa('Hash::Ordered')) {
      my @sps = ();
      for ( $qsp_files->keys()) {
         push @sps, short_species_name($_);
      }
      $self->set('querytaxa_short',\@sps);
   } else {
      die "query taxon inputpath is not a Hash::Ordered.\n";
   }
   for my $qsp ($qsp_files->keys()) {
      my $file = $qsp_files->get($qsp);
      my @gfiles = glob($file); # expand ~, *, ?, etc.
      die "query inputpath must specify unique file. $file expands to ", 
        scalar @gfiles, " files: ", join(" ", @gfiles), "\n" if(scalar @gfiles > 1);
      $file = abs_path($gfiles[0]); # expand ./ etc.
      $qsp_files->set($qsp, $file);
   }
   #  print join("\n", $qsp_files->values()), "\n";
   my $sp_files = $self->get('taxon_inputpath') // die "taxon_inputpath shash undefined.\n";
   for my $sp ($sp_files->keys()) {
      my $file = $sp_files->get($sp);
      my @gfiles = glob($file); # expand ~, *, ?, etc.
      die "query inputpath must specify unique file. $file expands to ", 
        scalar @gfiles, " files: ", join(" ", @gfiles), "\n" if(scalar @gfiles > 1);
      $file = abs_path($gfiles[0]); # expand ./ etc.
      $sp_files->set($sp, $file);
   }
   #  print join("\n", $sp_files->values()), "\n"; 

   my ($date, $time) = date_time();
   my $qspecies = $self->get('querytaxa_short');
   my $n_species = $self->get('n_species');
   my $filename_head = join('', @$qspecies) . "_vs_" . "$n_species" . "_" . $date;
   #  print STDERR "filename head: $filename_head \n";
   my $progress_filename = $filename_head . '.progress';
   #  print STDERR "progress_filename: $progress_filename \n";
   open my $fh_progress, ">", "$progress_filename";
   $self->set('fh_progress', $fh_progress);
   my $all_species_fasta_filename = "all_" . $n_species . "_species_" . $date . ".fasta";
   $self->set('all_species_fasta_filename', $all_species_fasta_filename);

   

}

1;


