package Overlap;
use strict;
use List::Util qw ( min max sum );

my $NO_ID = 'NOT_AN_ACTUAL_ID';

sub  new {
  my $class = shift;
  my $arg = shift; # either filename or string with contents of fasta file
  my $fraction = shift || 0.8;
  my $seed = shift || undef;
  my $args= {};
  my $self = bless $args, $class;

  if (defined $seed) {
    srand($seed);
  } else {
    srand();
  }
  my @ids = ();
  my %id_overlapseq = ();
  my %id_sequence = ();
  my @lines = ();
  if (! $arg =~ /\n/ and -f $arg) { # $arg is filename
    open my $fhin, "<$arg";
    @lines = <$fhin>;
    close $fhin;
  } else {			# treat arg as string
    @lines = split("\n", $arg);
  }
  my $id = $NO_ID;
  my $sequence = '';
  while (@lines) {
    my $line = shift @lines;
    if ($line =~ /^>/) {
      if ($id ne $NO_ID) {
	$id_sequence{$id} = $sequence;
	push @ids, $id;
      }
      $id = $line;
      $id =~ s/^>\s*//;
      $id =~ s/\s+$//;
      $sequence = '';
    } else {
      $line =~ s/^\s+//;
      $line =~ s/\s+$//;
      $sequence .= $line;
    }
  }
  if (! exists $id_sequence{$id}  and  $sequence ne '') { # take care of the last id-sequence pair.
    $id_sequence{$id} = $sequence;
    push @ids, $id;
  }
  $self->{id_seq} = \%id_sequence;
  $self->{ids} = \@ids;

  my $n_sequences = scalar @ids;
  my $seq_length = length $id_sequence{$ids[0]};
  $self->{align_length} = $seq_length;
  $self->{n_sequences} = $n_sequences;
  my @position_counts = ((0) x $seq_length);
  foreach my $id (@ids) {
    my $sequence = $id_sequence{$id};
    my $seql = length $sequence;
    die "Non-equal sequence lengths in alignment: $seq_length, $seql. id: $id \nsequence: $sequence\n" if($seql ne $seq_length); 
    for (my $i = 0; $i < $seq_length; $i++) {
      if (substr($sequence, $i, 1) ne '-') {
	$position_counts[$i]++;

      }
    }
  }
  $self->{position_counts} = \@position_counts;
  my $n_required = ($fraction >= 1)? $n_sequences: int ($fraction * $n_sequences) + 1;
  $self->{n_required} = $n_required;
  my $overlap_length = 0;
  my %id_overlapnongapcount = ();
  foreach my $position (0..@position_counts-1) {
    my $count = $position_counts[$position];
    if ($count >= $n_required) {
      $overlap_length++;
      foreach my $id (@ids) {
	my $char = substr($id_sequence{$id}, $position, 1);
	$id_overlapseq{$id} .= $char;
	$id_overlapnongapcount{$id}++ if($char ne '-');	
      }
    }
  }

  $self->{id_overlapseq} = \%id_overlapseq;
  $self->{id_overlapnongapcount} = \%id_overlapnongapcount;
  die "overlap length inconsistency??? $overlap_length \n" if($overlap_length != length $id_overlapseq{$ids[0]});
  $self->{overlap_length} = $overlap_length;
  #	$self->{ids} = \@ids;

  return $self;
}

sub align_fasta_string{
  my $self = shift;
  my $spacer = shift || '';
  my $align_fasta = '';
  foreach my $id (@{$self->{ids}}) {
    my $sequence = $self->{id_seq}->{$id};
    $align_fasta .= ">$spacer$id\n$sequence\n";
  }
  chomp $align_fasta;
  return $align_fasta;
}

sub overlap_fasta_string{
  my $self = shift;
  my $spacer = shift || '';
  my $overlap_fasta = '';
  foreach my $id (@{$self->{ids}}) {
    my $sequence = $self->{id_overlapseq}->{$id};
    $overlap_fasta .= ">$spacer$id\n$sequence\n";
  }
  chomp $overlap_fasta;
  return $overlap_fasta;
}

sub bootstrap_overlap_fasta_string{
  my $self = shift;
  my $spacer = shift || '';
  my %id_bootstrapoverlapseq = ();
  my $overlap_length = $self->{overlap_length};

  my @indices = ();
  for (1..$overlap_length) {
    my $index = int( rand($overlap_length) );
    push @indices, $index;
    #	$index_count{$index}++;
  }

  for my $id (@{$self->{ids}}) {
    my $std_overlap = $self->{id_overlapseq}->{$id};
    my $string = '';
    foreach my $index (@indices) {
      $string .=  substr($std_overlap, $index, 1);
    }
    $id_bootstrapoverlapseq{$id} = $string;
  }

  my $bofstring = '';
  foreach my $id (@{$self->{ids}}) {
    #	print "id, seq: $id; $sequence\n";
    my $sequence = $id_bootstrapoverlapseq{$id};
    $bofstring .= ">$spacer$id\n$sequence\n";
  }
  chomp $bofstring;
  return $bofstring;
}

sub get_overlap_length{
  my $self = shift;
  return $self->{overlap_length};
}
sub set_overlap_length{
  my $self = shift;
  $self->{overlap_length} = shift;
}


1;
