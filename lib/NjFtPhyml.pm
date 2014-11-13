package TomfyMisc;
use Exporter qw 'import';
@EXPORT_OK = qw 'run_quicktree';
use strict;


sub store_gg_info {	 #xx    # read in gene-genome association file
  # store in hash (keys: seqids, values: species). e.g. ('ATxxxx' => Arabidopsis_thaliana');
  my $gg_filename   = shift;
  my %seqid_species = ();
  if ( defined $gg_filename ) {
    if ( -f $gg_filename ) {
      open my $fh_gg, "<", "$gg_filename";
      while (<$fh_gg>) {
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
  return \%seqid_species;
}

sub run_quicktree {
  my $overlap_fasta_string       = shift;
  my $correction                 = shift || 'kimura';
  my $tmp_overlap_fasta_filename = "PID" . $$ . "_tmp_overlap_fasta";
  open my $fhtmp, ">", "$tmp_overlap_fasta_filename";
  print $fhtmp $overlap_fasta_string, "\n";
  close $fhtmp;

  my $tmp_overlap_stockholm_filename = "PID" . $$ . "_tmp_overlap_stockholm";
  system
    "sreformat stockholm $tmp_overlap_fasta_filename > $tmp_overlap_stockholm_filename";
  my $newick_out;
  if ( $correction eq 'kimura' ) {
    $newick_out = `quicktree -kimura $tmp_overlap_stockholm_filename`;
  } else {
    $newick_out = `quicktree  $tmp_overlap_stockholm_filename`;
  }

  $newick_out =~ s/\s+//g;	# remove whitespace
  $newick_out =~ s/;$//;
  return $newick_out;
}

sub run_fasttree {
  my $overlap_fasta_string = shift; # required
  my $fasttree_cl          = shift; # required 
  my $intree               = shift; # newick string for starting tree
  if ($intree) {
    open my $fh, ">", "tmp_intree";
    print $fh $intree, "\n";
    close $fh;
    $fasttree_cl .= " -intree tmp_intree ";
  }

  #print STDERR "fasttree cl: $fasttree_cl \n"; exit;
  my $fasttree_newick_out = "ft_newick_default_output";
  my $fasttree_stderr_out = "ft_stderr_default_output";
  my $fasttree_out = 'ft_default_output';

    my $temp_file = "PID" . $$ . "_ft_input_tmpfile";
    open my $fh, ">", "$temp_file";
    print $fh $overlap_fasta_string, "\n"; close $fh;
    $fasttree_out = `$fasttree_cl $temp_file 2>&1`; #  $stderr_outfile`;

    my @output_lines = split("\n", $fasttree_out);
    $fasttree_newick_out = pop @output_lines;
    my $ft_cpu_time = 
      ($output_lines[scalar @output_lines-1] =~ /Total \s+ time: \s+ ([0-9.]+) \s+ seconds/x)? $1 : '---';
    my $fasttree_loglikelihood = 
      ($output_lines[scalar @output_lines-2] =~ /Gamma [(] \d+ [)] \s+ LogLk \s+ = \s+ (-[0-9.]+) \s/x)? $1 : '---';

    my $stderr_outfile  = "PID" . $$ . "ft.stderr";
    open my $fh_stderr, ">", "$stderr_outfile" or die "couldnt open $stderr_outfile for writing.\n";
    print $fh_stderr join("\n", @output_lines), "\n";
    close $fh_stderr;

    # Gamma(20) LogLk = -5372.474 alpha = 1.562 rescaling lengths by 1.044   # parse ll out of ft stderr output.
    # my $fasttree_loglikelihood =
    #   ( $fasttree_stderr_out =~
    #     / Gamma [(] \d+ [)] \s+ LogLk \s+ = \s+ ([-] \d+ [.] \d*) \s+ alpha/xm
    #   ) ? $1 : undef;
    $fasttree_newick_out =~ s/\s+//g;
    $fasttree_newick_out =~ s/;$//;
    return ( $fasttree_newick_out, $fasttree_loglikelihood, $ft_cpu_time);
}


sub run_phyml{
  my $overlap_fasta_string = shift;
  my $initial_tree_newick = shift;
my $optimize = shift || 'tlr';
  my $phyml_obj = Phyml->new({
				#		    'dna_or_protein' => $dna_or_protein, # 'dna' or 'protein' (default is protein)
			      'fasta_string' => $overlap_fasta_string,
			      'optimize_param' => $optimize,
			      'n_bootstrap' => 0,
			      'initial_tree_newick' => $initial_tree_newick,
			      #			      'initial_tree_newick_file' => $initial_tree_newick_filename,
				#		    'initial_tree_newick_file' => $initial_tree_newick_file,
				#		    'alpha' => $alpha,
				#		    'p_invariant' => $p_invariant,
			      'subst_model' => 'WAG',
			      'n_rate_classes' => 4, 
			      'p_invariant' => 'e', # estimate the fraction of invariant sites
			     });

  $phyml_obj->run();

  #print "# phyml command line: ", $phyml_obj->{phyml_command_line}, "\n";
  #print $phyml_obj->{phyml_stdout}, "\n";

  my $phyml_loglikelihood = $phyml_obj->{log_likelihood};
  my $phyml_newick_out = $phyml_obj->{newick_out};
  my $phyml_cpu_time = $phyml_obj->{cpu_time};
  return ( $phyml_obj->{phyml_command_line}, $phyml_newick_out, $phyml_loglikelihood, $phyml_cpu_time );
}

1;
