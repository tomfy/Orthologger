package PipeSegGG;
use strict;

use parent 'PipelineSegment';

use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use Cwd qw(getcwd abs_path);

use TomfyMisc qw (date_time split_fasta_file split_fasta_string fix_fasta_files);

sub new {
   my $class = shift;
   my $segment_name = shift;
   my $pl_obj = shift;
   print STDERR "top of PipeSegGG constructor.\n";
   my $self  = $class->SUPER::new(); # bless $args, $class;
   $self->init(
               'segment_name' => $segment_name,
               'pipeline' => $pl_obj,
              );
   $self->hide('gg_string', 'seqid_species'); # add these to the set of hidden attributes.
   print STDERR "bottom of PipeSegGG constructor.\n"; #exit;
   return $self;
}

sub run{
   my $self = shift;
   my $taxon_filename = shift; # either an ordered hash with taxon:filename pairs, or a gg filename.

   my $pl_obj = $self->get('pipeline');
   $taxon_filename = $pl_obj->get('taxon_inputpath') if(!defined $taxon_filename); # default: get it from pipeline object.

   my ($date, $time) = date_time();
   my $gg_path =  abs_path($pl_obj->get('n_species') . "species_" . $date . ".gg");
   my ($v, $out_dir, $gg_filename) = File::Spec->splitpath($gg_path);
   $self->set(output_dir => $out_dir);
   $self->set(gg_filename => $gg_path);
   $pl_obj->set(gg_filename => $gg_path);

   if (blessed $taxon_filename and $taxon_filename->isa('Hash::Ordered')) { # taxon/fasta_file ordered hash
      if (scalar $taxon_filename->keys() > 0) {
         $self->construct_gg_from_fasta_files($taxon_filename);
      } else {
         die "in PipeSeqGG constructor; no taxon-fasta_file pairs specified.\n";
      }
   } else {
      die "PipeSegGG constructor failed. Optional argument should be taxon/file Hash::Ordered.\n";
   }
   $pl_obj->set(seqid_species => $self->get('seqid_species'));

   open my $fh_gg, ">", "$gg_path" or die "couldn't open $gg_path for writing.\n";
   print $fh_gg $self->get('gg_string');
   close $fh_gg;

   $self->set(state => 'complete');
 #  $self->set_successor_state('next');
   print $self->summary_string();
}

sub awaken{
   my $self = shift;
   my $id_species = $self->get('seqid_species');
   if (defined $id_species  and  blessed $id_species  and  $id_species->isa('Hash::Ordered')) {
      # don't need to do anything.
      print STDERR "gg awaken. do nothing branch.\n";
   } else {
      my $gg_filename = $self->get('gg_filename');
      if (-f $gg_filename) {
         $self->construct_gg_from_gg_file($gg_filename);
         $self->get('pipeline')->set(seqid_species => $self->get('seqid_species'));
         print STDERR "gg awaken. do something branch.\n";
      } else {
         die "Couldn't awaken gg segment, because gg file $gg_filename not present.\n";
      }
   }
   return;
}

sub construct_gg_from_fasta_files{
   my $self = shift;
   my $taxon_file = shift; # Hash::Ordered obj containing taxon/fasta_filename hash
   my %seqid_species = ();
   my $gg_string = '';
   for my $species ($taxon_file->keys()) {
      my $file = $taxon_file->get($species);
      $gg_string .= "$species:";
      my @ex_files = glob $file;
      $file = $ex_files[0];
      open my $fh, "<", "$file" or die "In construct_gg_from_fasta_files, couldn't open file [$file] for reading.\n";
      while (<$fh>) {
         if (/^>\s*(\S+)/) {    # this is a line with a sequence id.
            $gg_string .= " $1";
            $seqid_species{$1} = $species;
         }
      }
      $gg_string .= "\n";
      close $fh;
   }
   $self->set('seqid_species' => \%seqid_species);
   $self->set('gg_string' => $gg_string);
   return;
}

sub construct_gg_from_gg_file { # read in gene-genome association file
   # store in hash (keys: seqids, values: species). e.g. ('ATxxxx' => Arabidopsis_thaliana');
   my $self = shift;
   my $gg_filename   = shift;
   my %seqid_species = ();
   my $gg_string = '';
   if ( defined $gg_filename ) {
      if ( -f $gg_filename ) {
         open my $fh_gg, "<", "$gg_filename";
         while (<$fh_gg>) {
            $gg_string .= $_;
            my @cols = split( " ", $_ );
            my $species = shift @cols;
            $species =~ s/:$//; # remove final colon if present.
            for (@cols) {
               if ( exists $seqid_species{$_} ) {
                  warn "key $_ already stored with species: ",
                    $seqid_species{$_}, "\n";
               } else {
                  $seqid_species{$_} = $species;
               }
            }
         }
         close $fh_gg;
      } else {      # done storing gg_file info in hash %seqid_species
         die "$gg_filename: no such file.\n";
      }
   } else {
      die "gg filename undefined. \n";
   }
   $self->set('seqid_species' => \%seqid_species);
   $self->set('gg_string' => $gg_string);
   return;
}

sub summary_string{
   my $self = shift;
   my $summary_string = '';
   my @gg_lines = split("\n", $self->get('gg_string'));
   my $total = 0;
   for (@gg_lines) {
      my @xs = split(" ", $_);
      my $n_sequences = scalar @xs - 1;
      $summary_string .= sprintf("%30s  %8i  sequences\n", $xs[0], $n_sequences);
      $total += $n_sequences;
   }
   $summary_string .= sprintf("%30s  %8i  sequences\n", 'total', $total);
   $summary_string .= sprintf("gg file: %s\n", $self->get('gg_filename'));
   return $summary_string;
}

1;
