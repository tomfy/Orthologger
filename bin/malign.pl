#!/usr/bin/perl -w
use strict;

# read in a file with multiple families of sequences
# The format is for each family:
# 1) 1 line of info about family as whole
# 2) fasta for all sequences in family
# 3) one blank line
# for example: 
# Id Medtr1g004960.1 family. fam_size: 5 Medicago_truncatula
# >Medtr1g004960.1 
# MENSCMISAASRNLWKPAFQGSVNKTPRQIEREGVRNASNNNDDAKTNFCDKRKEKNIQTNEESDSSSSLEDEELRLHYKIHRSFNIPTIWGRHLKENSVPPKHNDCGDVTHDVFMHEIADILAYKLPYSVEHEKGLKGQSV*
# >Medtr3g051375.1 
# MKELDLDETTCHPNEVVVQKCMTSKKSRRKGRKTQVPWDERCKNVEYSPKNKRRTKKGMWNVSKNNDDVEFNALHKGKEKMIQTNVESNSSSFSEDEELRLALQKIVDLSISTKIWGRHLKVFARRN*
# >Medtr1g026650.1 
# MPGNNLMALLMMIFISSSTMSCIWPFISNVDLLSSSVISLIKQINKHLIVESICGGSGPIGGEVSLKPRNPQVPWDERCKNAENLVIKTKEEPKVSDIESNLVGGIVISEPRKIEIEGMWNVTNNNDQVESNFRHKGKEEMIQTNEESDSSFLHSLLKMRS*
# >Medtr3g453170.1 
# MKELDLDETIYHPNEVVVQKCMTSKKSRRKGRKTQVAWDERCKNAESSPKNKRRTRKGMWNVSKNNDDDEFNALHKGKEKMIQTNVESNSSSFSEDEELRLALQNSLIFQYPLKFGEDTSKSLLEGTNSRHQFNIEEGIKENSFEYFHIPNNIRINTLNL*
#
# Id Medtr...
my $input_filename = shift;
my $maxiters = shift || 2; # 16;
# my $n_taxa = shift || 21;
# my $min_taxa = 4;
my $state = 'idline';		# other possible values: fasta
my ($qid, $fam_size, $taxa, $idline, $fasta, $do);
open my $fh_in, "<", "$input_filename";
my $tmp_filename = $input_filename . "_tmp";
my $stderr_filename = $input_filename . ".stderr";
while (<$fh_in>) {
  if ($state eq 'idline') {
    my @cols = split(" ", $_);
    ($qid, $fam_size, $taxa) = @cols[1,4,5];
    my @species = split(",", $taxa);
    $do = 1;
    # if(scalar @species < $min_taxa){ # need to have at least 4 taxa (Medtr + 3 monocots)
    #   $do = 0;
    # }elsif(scalar @species >= $n_taxa-1){ # if have at least 20 out of 21, guarantees 3 monocots
    #   $do = 1;
    # }else{ # need to look at taxa in more detail:
    #   my ($monocot_count, $selaginella_present) = check_taxon_list($taxa);
    #   $do = ($monocot_count >= 3  and  $Selaginella_present)? 1 : 0;
    # }
    $idline = $_;
    $state = 'fasta';
    $fasta = '';
  } elsif ($state eq 'fasta') {
    if(/^>/){
      $fasta .= "\n" . $_;
    }else{
    my $fasta_line = $_;
    chomp $fasta_line;
    $fasta .= $fasta_line;
#    print "fasta: $fasta \n";
#exit;
  }
    if (/^\s*$/) {
      chomp $idline;
      print "$idline  $do \n";
      if($do){
      open my $fh, ">", "$tmp_filename";
      $fasta =~ s/^\n//; # not global - just initial one if present.
      print $fh "$fasta \n";
      close $fh;
      my $fasta_alignment  = 
	`muscle -in $tmp_filename -maxiters $maxiters 2> $stderr_filename`;
#    `mafft --retree 2 --inputorder $tmp_filename 2> $stderr_filename`; # mafft fastest
#  `mafft --auto  --inputorder $tmp_filename 2> $stderr_filename`; # mafft slower, better (presumably)
      print "$fasta_alignment \n";
    }else{
      print "\n";
    }
      $state = 'idline';
    }
  }
}


sub check_taxon_list{ # check list of taxa to see how many monocot species are
# present, and whether selaginella is present.
my $taxa = shift;

my $monocot_count = 0;
$monocot_count += 1 if($taxa =~ /Oryza_sativa/);
$monocot_count += 1 if($taxa =~ /Sorghum_bicolor/);
$monocot_count += 1 if($taxa =~ /Brachypodium_distachyon/);
$monocot_count += 1 if($taxa =~ /Zea_mays/);
return ($monocot_count, $taxa =~ (/Selaginella/)? 1 : 0);
}
