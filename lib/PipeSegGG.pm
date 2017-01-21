package PipeSegGG;
use strict;
use base 'PipelineSegment';
use List::Util qw (min max sum);

sub new {
   my $class = shift;
   my $gg_href_or_file = shift; # 
   my $args  = {};
   my $self  = bless $args, $class;

   if (ref $gg_href_or_file eq 'HASH' ) { # taxon/fasta_file hash
      if (scalar keys %$gg_href_or_file > 0) {
         $self->construct_gg_from_fasta_files($gg_href_or_file);
      } else {
         die "in PipeSeqGG constructor; no taxon-fasta_file pairs specified.\n";
      }
   } elsif (-f $gg_href_or_file) { # gg file
      $self->construct_gg_from_gg_file($gg_href_or_file);
   } else {
      die "PipeSegGG constructor failed. Argument should be taxon/file hashref, or gg file name.\n";
   }
   return $self;
}

sub construct_gg_from_fasta_files{
   my $self = shift;
   my $taxon_file = shift; # hashref; e.g. {'Arabidopsis_lyrata' => 'Alyrata.fasta'}
   my %seqid_species = ();
   my $gg_string = '';
   while ( my ($species, $file) = each %$taxon_file) {
      $gg_string .= "$species:";
      my @ex_files = glob $file;
      $file = $ex_files[0];
      open my $fh, "<", "$file" or die "In construct_gg_from_fasta_files, couldnt open file [$file] for reading.\n";
      while (<$fh>) {
         if (/^>\s*(\S+)/) {	# this is a line with a sequence id.
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
            $species =~ s/:$//;	# remove final colon if present.
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
      } else {	    # done storing gg_file info in hash %seqid_species
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
