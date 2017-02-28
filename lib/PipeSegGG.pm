package PipeSegGG;
use strict;

use base 'PipelineSegment';

use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use Cwd qw(getcwd abs_path);

use TomfyMisc qw (date_time split_fasta_file split_fasta_string fix_fasta_files);

sub new {
   my $class = shift;
   my $pl_obj = shift;
   print STDERR "top of PipeSegGG constructor.\n";
   my $self  = $class->SUPER::new($pl_obj); # bless $args, $class;
   print STDERR "bottom of PipSegGG constructor.\n";
   return $self;
}

sub run{
   my $self = shift;
   my $gg_ordhash_or_file = shift; # 

   my ($date, $time) = date_time();
   my $pl_obj = $self->get('pipeline');
   print STDERR "in PSGG. pl_obj defined? ", (defined $pl_obj)? 'yes' : 'no', "  $date ", $pl_obj->get('n_species'), "\n";
   my $gg_filename =  abs_path($pl_obj->get('n_species') . "species_" . $date . ".gg");
   
   $pl_obj->set('gg_filename',$gg_filename);

   $gg_ordhash_or_file = $pl_obj->get('taxon_inputpath') if(!defined $gg_ordhash_or_file); # i.e. default is get it from pipeline object.

   if (blessed $gg_ordhash_or_file and $gg_ordhash_or_file->isa('Hash::Ordered')){ # taxon/fasta_file ordered hash
      if (scalar $gg_ordhash_or_file->keys() > 0) {
         $self->construct_gg_from_fasta_files($gg_ordhash_or_file);
      } else {
         die "in PipeSeqGG constructor; no taxon-fasta_file pairs specified.\n";
      }
   } elsif (-f $gg_ordhash_or_file) { # gg file
      $self->construct_gg_from_gg_file($gg_ordhash_or_file);
   } else {
      die "PipeSegGG constructor failed. Argument should be taxon/file Hash::Ordered, or gg file name.\n";
   }
       $pl_obj->set('seqid_species', $self->get('seqid_species'));
       $self->get('hidden')->{gg_string} = 1; # hide gg_string, so don't include gg_string in stringified obj.

       open my $fh_gg, ">", "$gg_filename" or die "couldn't open $gg_filename for writing.\n";
       print $fh_gg $self->get('gg_string');
       close $fh_gg;
    }

     sub construct_gg_from_fasta_files{
        my $self = shift;
        my $taxon_file = shift; # Hash::Ordered obj containing taxon/fasta_filename hash
        my %seqid_species = ();
        my $gg_string = '';
        #   while ( my ($species, $file) = each %$taxon_file) {
        for my $species ($taxon_file->keys()) {
           my $file = $taxon_file->get($species);
           $gg_string .= "$species:";
           my @ex_files = glob $file;
           $file = $ex_files[0];
           open my $fh, "<", "$file" or die "In construct_gg_from_fasta_files, couldn't open file [$file] for reading.\n";
           while (<$fh>) {
              if (/^>\s*(\S+)/) { # this is a line with a sequence id.
                 $gg_string .= " $1";
                 $seqid_species{$1} = $species;
              }
           }
           $gg_string .= "\n";
           close $fh;
        }
        $self->set('seqid_species', \%seqid_species);
        $self->set('gg_string', $gg_string);
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
         } else {   # done storing gg_file info in hash %seqid_species
            die "$gg_filename: no such file.\n";
         }
      } else {
         die "gg filename undefined. \n";
      }
      $self->set('seqid_species', \%seqid_species);
      $self->set('gg_string', $gg_string);
      return;
   }

   sub summary{
      my $self = shift;
      my $summary_string = '';
      my @gg_lines = split("\n", $self->get('gg_string'));
      my $total = 0;
      for (@gg_lines) {
         #   print $_, "\n";
         my @xs = split(" ", $_);
         my $n_sequences = scalar @xs - 1;
         #     $summary_string .= $xs[0] . "  " . (scalar @xs - 1) . "  sequences.\n";
         $summary_string .= sprintf("%30s  %8i  sequences\n", $xs[0], $n_sequences);
         $total += $n_sequences;
      }
      $summary_string .= sprintf("%30s  %8i  sequences\n", 'total', $total);
      return $summary_string;
   }


   1;
