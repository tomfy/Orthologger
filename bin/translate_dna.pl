#!/usr/bin/perl -w
use strict;

#my $input_file = shift;

my %rnacodon_aa = (             # X: STOP
                   UUU => 'F', UUC => 'F',
                   UUA => 'L', UUG => 'L', CUU => 'L', CUC => 'L', CUA => 'L', CUG => 'L',
                   AUU =>'I', AUC => 'I', AUA => 'I',
                   AUG => 'M',
                   GUU => 'V', GUC => 'V', GUA => 'V', GUG => 'V',
                   UCU => 'S', UCC => 'S', UCA => 'S', UCG => 'S',
                   CCU => 'P', CCC => 'P', CCA => 'P', CCG => 'P',
                   ACU => 'T', ACC => 'T', ACA => 'T', ACG => 'T',
                   GCU => 'A', GCC => 'A', GCA => 'A', GCG => 'A',
                   UAU => 'Y', UAC => 'Y',
                   UAA => 'X', UAG => 'X', UGA => 'X',
                   CAU => 'H', CAC => 'H',
                   CAA => 'Q', CAG => 'Q',
                   AAU => 'N', AAC => 'N',
                   AAA => 'K', AAG => 'K',
                   GAU => 'D', GAC => 'D',
                   GAA => 'E', GAG => 'E',
                   UGU => 'C', UGC => 'C',
                   UGG => 'W',
                   CGU => 'R', CGC => 'R', CGA => 'R', CGG => 'R',
                   AGU => 'S', AGC => 'S',
                   AGA => 'R',
                   AGG => 'R',
                   GGU => 'G', GGC => 'G', GGA => 'G', GGG => 'G',
                  );

my %nt_c = (A => 'T', T => 'A', G => 'C', C => 'G');

my %codon_aa = ();

while (my($rc, $aa) = each %rnacodon_aa) {
   my $dc = $rc;
   $dc =~ s/U/T/g;
   $codon_aa{$dc} = $aa;
}


my @sorted_codons = sort {$codon_aa{$a} cmp $codon_aa{$b}} keys %codon_aa;

for my $c (@sorted_codons) {
   my $aa = $codon_aa{$c};
}

my %id_longtransorf = ();
my $trunc_id = undef;
while (<>) {
   if (/^>(\S+)/) {
      my $id = $1;
      $trunc_id = $id;
      $trunc_id =~ s/_(i\d+)\s*$//;
   #   print "$trunc_id  $1  ";
   } else {
      s/^[^a-zA-Z]+//;
      s/[^a-zA-Z]+$//;
      my $maxorfl = -1;
      my $Lmax = -1;
      my $longest_orfseq = undef;
   #   my @orflengths = ();    # longest orfs for the 6 reading frames.
      for my $strand ('forward', 'reverse') {
         for my $offset (0..2) {
            my $aa_sequence =  translate_dna($_, $offset, $strand);
            my ($L, $orfl, $orfseq) =  aas2lorf($aa_sequence);
         #   push @orflengths, $orfl;
            if ($orfl > $maxorfl) {
               $maxorfl = $orfl;
               $longest_orfseq = $orfseq;
               $Lmax = $L;
            }
         }
      }
   #   print "orf lengths: ", join(" ", @orflengths[0..5]), "   ";
  #    @orflengths = sort {$b <=> $a} @orflengths;
   #   print "sorted orf lengths: ", join(" ", @orflengths[0..5]), "\n";
  #    print "$Lmax  $maxorfl  $longest_orfseq \n";
      if (!exists $id_longtransorf{$trunc_id}) {
         $id_longtransorf{$trunc_id} = $longest_orfseq;
      } elsif (length $longest_orfseq > length $id_longtransorf{$trunc_id}) {
         $id_longtransorf{$trunc_id} = $longest_orfseq;
      }
   }
}

while(my ($id, $seq) = each %id_longtransorf){
   print ">$id  ", length $seq, "\n", "$seq \n";
}

sub translate_dna{
   my $dna_seq = shift;
   my $offset = shift // 0;
   my $strand = shift // 'forward';

   my $pep_seq = '';

   $dna_seq = reverse_complement($dna_seq) if($strand eq 'reverse');
   substr($dna_seq, 0, $offset, '');
   while (length $dna_seq >= 3) {
      my $codon = substr($dna_seq, 0, 3, '');
      if (1 or $strand eq 'forward') {
         $pep_seq .= $codon_aa{$codon} // '?';
      } elsif ($strand eq 'reverse') {
         my $revcomp_codon = reverse_complement($codon);
         $pep_seq = ($codon_aa{$revcomp_codon} // '?') . $pep_seq;
      } else {
         die "strand is $strand; should be forward or reverse.\n";
      }
   }
   return $pep_seq;
}

sub reverse_complement{
   my $dna = shift;
   my $rc_dna = '';
   my @chars = split('', $dna);
   for (@chars) {
      $rc_dna = $nt_c{$_} . $rc_dna;
   }
   return $rc_dna;
}


sub aas2lorf{
   my $aaseq = shift;

   my $id = 'No_Id';
   my $state = 'O';
   my $orf = '';
   my $long_orf = '';
#   my $position = 0;
#   my $start_pos = 0;
#   my $long_orf_start_pos = 0;
#   my $long_orf_end_pos = 0;
   $aaseq =~ s/\s+$//g;         # remove all whitespace
   my @chars = split('', $aaseq);

   for my $c (@chars) {
      if ($c eq 'M'  and  $state eq 'C') {
         $state = 'O';
         $orf = 'M';
      } elsif ($c eq 'X'  and  $state eq 'O') {
         if(length $orf > length $long_orf){
            $long_orf = $orf;
       #     $long_orf_start_pos = $start_pos;
       #     $long_orf_end_pos = $position;
         }
         $orf = '';
         $state = 'C';
      } else {
         if ($state eq 'O') {
            $orf .= $c;
         }
      }
    #  $position++;
   }
      if(length $orf > length $long_orf){
      $long_orf = $orf;
  #    $long_orf_start_pos = $start_pos;
 #     $long_orf_end_pos = $position;
   }
   return (length $aaseq, length $long_orf, $long_orf);
}


