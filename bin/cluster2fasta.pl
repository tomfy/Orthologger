#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $fasta_filename = shift || undef;

my $id_seq = store_fasta_in_hash($fasta_filename);

while (<>) {			# read clusterfile, 1 cluster per line
  s/^([^:]+):\s*//;
  my $cluster_name = $1;
  my @ids = split(" ", $_);
  my $out_string = '';
  my $out_filename = "$cluster_name.fasta";
  for (@ids) {
    $out_string .=  ">$_\n";
    if (exists $id_seq->{$_}) {
      $out_string .=  $id_seq->{$_} . "\n";
    } else {
      warn "No sequence found for id: $_\n";
    }
  }
  print $out_string;
  # open my $fh, ">", "$out_filename";
  # print $fh $out_string;
  # close $fh;
}


sub store_fasta_in_hash{
  my $fasta_filename = shift;

  die "File $fasta_filename does not exist." if(! -f $fasta_filename);

  my %id_seq = ();
  open my $fh, "<", "$fasta_filename";
  my $id = undef;
  my $sequence = undef;
  while (<$fh>) {
    if (/^>/) {
      if (defined $id) {
	$id_seq{$id} = $sequence;
      }
      /^>([^ |]+)/;
      $id = $1;
      $sequence = '';
    } else {
      /^\s*([^\s]+)/;
      $sequence .= $1;
    }
  }
  return \%id_seq;
}
