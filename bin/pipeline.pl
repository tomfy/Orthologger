#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $max_eval = 1e-12;		# default.
my $max_fam_size = 140;
my $ggfilename = undef;
my $m8_filename = undef;
my $fasta_infile = undef;
my $n_pieces = 8;
# my $abc_filename = undef;
my $first_of_species_factor = 1; # has no effect if no ggfile.
# by setting this to some large number (e.g. 1e20) will keep the best match
# from each species, if it is within this factor of the e-value threshold
# so if it is 1, it has no effect.

# my $best_model_only = 0;

# Process long cl options
GetOptions(
	   'input_filename=s' => \$m8_filename,
	   'max_eval=s' => \$max_eval,
	   'fam_size_limit=i' => \$max_fam_size,
	   'ggfilename=s' => \$ggfilename,
	   'fasta_infile=s' => \$fasta_infile,
	   # 'only_best_model=s' => \$best_model_only
	   'n_pieces=i' => \$n_pieces,
	  );
my $abc_filename = '';
if (defined $ggfilename) {
  if (-f $ggfilename) {
    $abc_filename = `m8toabc_fams.pl -input $m8_filename -max_eval $max_eval -fam_size_limit $max_fam_size -ggfile $ggfilename`;
  } else {
    die "No file with path $ggfilename. Exiting. \n";
  }
} else {
  $abc_filename = `m8toabc_fams.pl -input $m8_filename -max_eval $max_eval -fam_size_limit $max_fam_size`;
}
$abc_filename =~ s/\s+$//; 
print STDERR "[$abc_filename]\n";

my $abc_part_filenames = `split_abc.pl $abc_filename $n_pieces`;
my @abc_parts = split(" ", $abc_part_filenames);

print "[" ,join("][", @abc_parts), "]\n";

print "PID: $$ \n";
my @fastas_filenames = ();
die "fasta input file undefined or file does not exist.\n" unless(defined $fasta_infile and -f $fasta_infile); 
for my $the_abc_part (@abc_parts) {
  my $output_fastas_filename = $the_abc_part;
  $output_fastas_filename =~ s/abc$/fastas/;
  print "seq+... output filename: $output_fastas_filename \n";
  system("seq+matches2fasta.pl -gg $ggfilename -abc_file $the_abc_part -fasta_infile $fasta_infile -output_filename $output_fastas_filename");
  push @fastas_filenames, $output_fastas_filename
}

# fork processes to do alignment, tree finding.
my $n_alignments_to_do = scalar @fastas_filenames;
for my $a_fastas_filename (@fastas_filenames) {
  my $malign_out_filename = $a_fastas_filename;
  $malign_out_filename =~ s/fastas$/alfastas/;
  my $pid = fork(); # returns 0 to child process, pid of child to parent process.
  if ($pid == 0) {  # child process
    my $malign_stdout = `malign.pl  -input $a_fastas_filename  -align muscle  -quality quick  -output $malign_out_filename`;
    print "malign.pl finished aligning $a_fastas_filename; output file: $malign_out_filename. \n";
    my $output_newick_filename = $malign_out_filename;
    $output_newick_filename =~ s/alfastas$/newicks/;
    my $nj_ft_bs_stdout = `nj_ft_bs.pl -gg $ggfilename -input $malign_out_filename -output $output_newick_filename`;
    exit(0);
  }
}
my $children_done = 0;
while (wait() != -1) {
  $children_done++; print "Number of alignments finished: $children_done out of $n_alignments_to_do. \n";
}
print "Done waiting \n";
