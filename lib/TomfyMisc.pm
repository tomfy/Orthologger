package TomfyMisc;
use Exporter qw 'import';
@EXPORT_OK = qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring file_format_is_abc_iie';
use strict;

sub timestring{
  my $s_time = shift; # time in seconds (e.g. from time() )
  my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  my @ltime = localtime($s_time);
  my $t_string = $ltime[2] . ":" . $ltime[1] . ":" . $ltime[0];
my $d_string = $ltime[5]+1900 . "-" . $months[$ltime[4]] . "-" . $ltime[3];
return "$t_string, $d_string.";
}

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
  #  print STDERR "fasttree cl: ", $fasttree_cl, "\n";
  #  print "AAAAA: ft newickout, stderrout: $fasttree_newick_out, $fasttree_stderr_out \n";
  if (0) {		# doesn't seem to work now, did before. Why???
    my $run3_return_value = run3(
				 "$fasttree_cl",        \$overlap_fasta_string,
				 \$fasttree_newick_out, \$fasttree_stderr_out
				);
    #    print STDERR "in run_fasttree, run3 return value: [", $run3_return_value, "]\n";
  } else {			# do using a temporary file
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

    #  print STDERR "FT OUT: \n[ ", $fasttree_newick_out, " ]\n\n"; 
    #  print STDERR "ft cputime:  $ft_cpu_time.  ft ln(L): $fasttree_loglikelihood \n";

    #print STDERR "FT stderr out: \n", $fasttree_newick_out, "\n";
    #print STDERR "FT OUT stderr: \n", $fasttree_stderr_out, "\n\n";

    # Gamma(20) LogLk = -5372.474 alpha = 1.562 rescaling lengths by 1.044   # parse ll out of ft stderr output.
    # my $fasttree_loglikelihood =
    #   ( $fasttree_stderr_out =~
    #     / Gamma [(] \d+ [)] \s+ LogLk \s+ = \s+ ([-] \d+ [.] \d*) \s+ alpha/xm
    #   ) ? $1 : undef;
    $fasttree_newick_out =~ s/\s+//g;
    $fasttree_newick_out =~ s/;$//;
    return ( $fasttree_newick_out, $fasttree_loglikelihood, $ft_cpu_time);
  }
}

sub run_phyml{
  my $overlap_fasta_string = shift;
  my $initial_tree_newick = shift;
  my $to_optimize = shift || 'tlr'; # or 'lr'
  my $phyml_obj = Phyml->new({
				#		    'dna_or_protein' => $dna_or_protein, # 'dna' or 'protein' (default is protein)
			      'fasta_string' => $overlap_fasta_string,
			      'optimize_param' => $to_optimize,
			      'n_bootstrap' => 0,
			      'initial_tree_newick' => $initial_tree_newick,
			      #			      'initial_tree_newick_file' => $initial_tree_newick_filename,
				#		    'initial_tree_newick_file' => $initial_tree_newick_file,
				#		    'alpha' => $alpha,
				#		    'p_invariant' => $p_invariant,
			      'subst_model' => 'WAG',
			      'n_rate_classes' => 4, 
			      'p_invariant' => 'e',
			     });

  $phyml_obj->run();

  #print "# phyml command line: ", $phyml_obj->{phyml_command_line}, "\n";
  #print $phyml_obj->{phyml_stdout}, "\n";

  my $phyml_loglikelihood = $phyml_obj->{log_likelihood};
  my $phyml_newick_out = $phyml_obj->{newick_out};
  my $phyml_cpu_time = $phyml_obj->{cpu_time};
  return ( $phyml_obj->{phyml_command_line}, $phyml_newick_out, $phyml_loglikelihood, $phyml_cpu_time );
}

sub file_format_is_abc{
   my $filename = shift;
   my $ok_line_count = 0;
   open my $fh, "<", "$filename" or die "couldn't open $filename for reading.\n";
   for (1..3) {
      my $line = <$fh>;
      last unless($line =~ /^\S+\s+\S+\s+\S+/);
      $ok_line_count++;
   }
   close $fh;
   return ($ok_line_count == 3);
}

sub file_format_is_iie{
   my $filename = shift;
   my $ok_line_count = 0;
   open my $fh, "<", "$filename" or die "couldn't open $filename for reading.\n";
   my $line = <$fh>;
   $ok_line_count++ if($line =~ /^\S+/);
   $line = <$fh>;
   last unless($line =~ /^\s+\S+\s+\S+/);
   $ok_line_count++;
   close $fh;
   return ($ok_line_count == 2);
}

sub file_format_is_abc_iie{
   my $filename = shift;
   my $format = 'other';
   if (file_format_is_abc($filename)) {
      $format = 'abc';
   } elsif (file_format_is_iie($filename)) {
      $format = 'iie';
   }
   return $format;
}

1;
