#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# read in a file with multiple families of sequences
# The format ("fastas" format) is for each family:
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
my $start_time_seconds = time();
my $input_filename;
my $align_program = 'muscle';
my $quality = 'medium';
open my $fh_out, '>&', STDOUT or die "$!";
my $output_file = undef;
my $maxiters;

# Process long cl options
GetOptions(
	   'input_filename=s'      => \$input_filename,
	   'align_program=s'      => \$align_program, # 'muscle' or 'mafft'
	   'quality=s' => \$quality, # 'best' or 'quick'
	   'output=s' => \$output_file # default is to write output to STDOUT
	  );
if (defined $output_file) {
  close $fh_out;
  open $fh_out, ">", "$output_file" or die "$!";
}
# my $n_taxa = shift || 21;
# my $min_taxa = 4;
#print STDERR "[$align_program]  [$quality]\n";
print STDERR "[$align_program] [$input_filename] [$quality]\n"; #exit;

my $state = 'idline';		# other possible values: fasta
my ($qid, $fam_size, $taxa, $idline, $fasta, $do);
open my $fh_in, "<", "$input_filename" or die "Couldnt open $input_filename for reading.\n";
my $tmp_filename = $input_filename . "_tmp";
my $stderr_filename = $input_filename . ".stderr";
my $count_alignments_done = 0;

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
    if (/^>/) {
      $fasta .= "\n" . $_;
    } else {
      my $fasta_line = $_;
      chomp $fasta_line;
      $fasta .= $fasta_line;
      #    print "fasta: $fasta \n";
      #exit;
    }
    if (/^\s*$/) {		# only whitespace
      chomp $idline;
      print $fh_out "$idline  $do \n";
      if ($do) {
	open my $fh_temp, ">", "$tmp_filename";
	$fasta =~ s/^\n//; # not global - just initial one if present.
	print $fh_temp "$fasta \n";
	close $fh_temp;

	my $alignment_cl;
	if ($quality ne 'best'  and  $quality  ne  'quick' and $quality ne 'medium') {
	  warn "Parameter quality set to $quality. Not implemented. Using default ('medium') quality.";
	}
	if ($align_program eq 'muscle') {
	  if ($quality eq 'best') {
	    $maxiters = 25;
	}elsif($quality eq 'medium'){
	$maxiters = 4;
	  } else {		# ($quality eq 'quick'){
	    $maxiters = 2;
	  }
	  $alignment_cl = "muscle -in $tmp_filename -maxiters $maxiters -maxhours 0.4  2> $stderr_filename";
	} elsif ($align_program eq 'mafft') {
	  if ($quality eq 'best') {
	    $alignment_cl = 
              "mafft --maxiterate 25 --localpair  --inputorder $tmp_filename 2> $stderr_filename"; # L-INS-i (best mafft);
         } elsif ($quality eq 'groovy'){ # better than 'medium', not as good as 'best'
             "mafft --retree 2 --maxiterate 25 --inputorder $tmp_filename 2> $stderr_filename";
} elsif ($quality eq 'medium'){
$alignment_cl = "mafft --retree 2 --maxiterate 2 --inputorder $tmp_filename 2> $stderr_filename"; # "mafft --auto $tmp_filename 2> $stderr_filename";  
} else {		# ($quality eq 'quick'){
	    $alignment_cl = "mafft --retree 2 --inputorder $tmp_filename 2> $stderr_filename"; # mafft quick
	  }
	} else {
	  warn "Parameter align_program set to $align_program. Not implemented. Using default: 'muscle'."
	}
	print STDERR "Alignment command line: [$alignment_cl] \n";
	my $fasta_alignment = `$alignment_cl`; 
	$count_alignments_done++;
	print $fh_out "$fasta_alignment \n";
      } else {
	print $fh_out "\n";
      }
      $state = 'idline';
    }
  }
}
my $end_time_seconds = time();
print STDERR "Alignments done: $count_alignments_done in ", $end_time_seconds - $start_time_seconds, " seconds\n";


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
