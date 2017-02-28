package PhylogenomicPipeline;
use strict;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use File::Basename 'dirname';
#use File::Spec qw(splitpath);
use Cwd qw(abs_path getcwd);

use parent 'Pipeline';

my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir); # collapses the bin/../lib to just lib
}
use lib $libdir;

use TomfyMisc qw(date_time short_species_name);
use PipeSegGG;
use PipeSegBlast;
use PipeSegRbmcl;
use PipeSegFamilies;

sub new {
   my $class = shift;
   my $control_filename = shift; #
   my $args = shift;
   my $default_params = Hash::Ordered->new(
                                           attributes => [],
                                           hidden => {},
                                           min_sequence_length => 20,
                                           n_pieces => 3,
                                          );
   print STDERR "SSS: ", join(" ", @$args), "\n";
   $default_params->merge(@$args);
   # for(@$args){ # $args override the default parameters.
   #    $default_params->{$_} = $args->{$_};
   # }
   my $self = $class->SUPER::new();
   
   print STDERR "In PhylogenomicPipeline constructor\n";
   #  my $self  = bless $args, $class;
   $self->initialize($control_filename, $default_params);         # $params);
   return $self;
}

sub initialize{
   my $self = shift;

     my $control_filename = shift;
      my $params = shift;
      $self->SUPER::initialize($control_filename, $params);
print STDERR "Top of PhylogenomicPipeline::initialize() \n";

   $self->set('n_species',  scalar $self->get('taxon_inputpath')->keys()); # scalar keys %{$self->get('taxon_inputpath')});

   my $qsp_files = $self->get('querytaxon_inputpath') // die "querytaxon_inputpath Hash::Ordered undefined.\n";
   if (blessed $qsp_files and $qsp_files->isa('Hash::Ordered')) {
      #   $self->set('querytaxa', $qsp_files->keys());
      my @sps = ();
      for ( $qsp_files->keys()) {
         #  for (@sps) {
         #         print STDERR "$_  ", short_species_name($_), "\n";
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
   print join("\n", $qsp_files->values()), "\n";
   # exit;
   my $sp_files = $self->get('taxon_inputpath') // die "taxon_inputpath shash undefined.\n";
   #   while (my ($sp, $file) = each %$sp_files) {
   for my $sp ($sp_files->keys()) {
      my $file = $sp_files->get($sp);
      my @gfiles = glob($file); # expand ~, *, ?, etc.
      die "query inputpath must specify unique file. $file expands to ", 
        scalar @gfiles, " files: ", join(" ", @gfiles), "\n" if(scalar @gfiles > 1);
      $file = abs_path($gfiles[0]); # expand ./ etc.
      $sp_files->set($sp, $file);
   } 
   print join("\n", $sp_files->values()), "\n"; 
   # exit;

   my ($date, $time) = date_time();
   my $qspecies = $self->get('querytaxa_short');
   my $n_species = $self->get('n_species');
   my $filename_head = join('', @$qspecies) . "_vs_" . "$n_species" . "_" . $date;
   print STDERR "filename head: $filename_head \n";
   my $progress_filename = $filename_head . '.progress';
   print STDERR "progress_filename: $progress_filename \n";
   open my $fh_progress, ">", "$progress_filename";
   $self->set('fh_progress', $fh_progress);
   my $all_species_fasta_filename = "all_" . $n_species . "_species_" . $date . ".fasta";
   $self->set('all_species_fasta_filename', $all_species_fasta_filename);

}


1;


