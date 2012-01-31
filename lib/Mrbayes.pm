package Mrbayes;
use strict;
use List::Util qw ( min max sum );
use TlyUtil qw ( Kolmogorov_Smirnov_D );

# this is an object to facilitate running MrBayes (Bayesian phylogeny
# inference program). The main functionality it adds is a relatively
# easy way to control the criteria for deciding when to end a run.

my $default_chunk_size = 2000;

sub  new {
  my $class = shift;
  my $arg = shift;	       # a hashref for setting various options
  my $default_arguments = {'alignment_nex_filename' => undef,
			   'file_basename' => undef,
			   'seed' => undef,
			   'swapseed' => undef,
			   'n_runs' => 2,
			   'n_temperatures' => 4,
			   'temperature_gap' => 0.25,
			   'chunk_size' => $default_chunk_size,
			   'print_freq' => undef,
			   'sample_freq' => 20,
			   'burnin_frac' => 0.1,
			   'diag_freq' => undef,
			   'converged_chunks_required' => 10,

			   'fixed_pinvar' => undef, # undef -> leaves default of uniform(0,1) in effect
			   # convergence criteria
			   'splits_min_hits' => 25,
			   'splits_max_ok_stddev' => 0.03,
			   'splits_max_ok_avg_stddev' => 0.01,
			   'modelparam_min_ok_ESS' => 250,
			   'modelparam_max_ok_PSRF' => 1.02,
			   'modelparam_max_ok_KSD' => 0.2,
			   'ngens_run' => 0
			  };
  my $self = bless $default_arguments, $class;

  foreach my $option (keys %$arg) {
    warn "Unknown option: $option in Mrbayes constructor.\n" if(!exists $self->{$option});
    $self->{$option} = $arg->{$option};
  }
  $self->{print_freq} = $self->{chunk_size} if(!defined $self->{print_freq});
  $self->{diag_freq} = $self->{chunk_size} if(!defined $self->{diag_freq});
  # print "print, diag freq: ", $self->{print_freq}, "  ", $self->{diag_freq}, "\n";

  my $alignment_nex_filename = $self->{alignment_nex_filename};
  #  if(defined $self->{file_basename}){
  #  my  $file_basename;
  if (!defined $self->{file_basename}) {
    my $file_basename = $alignment_nex_filename;
    $file_basename =~ s/[.]nex$//; # delete .nex ending
    $self->{file_basename} = $file_basename;
  }
  my $n_runs = $self->{n_runs};
  my $burnin_frac = $self->{burnin_frac};
  my $n_temperatures = $self->{n_temperatures};
  my $temperature_gap = $self->{temperature_gap};
  my $sample_freq = $self->{sample_freq};
  my $print_freq = $self->{print_freq};
  my $chunk_size = $self->{chunk_size};
  my $fixed_pinvar = $self->{fixed_pinvar};
  my $prset_pinvarpr = (defined $fixed_pinvar)? "prset pinvarpr=fixed($fixed_pinvar);\n" : '';

  my $begin_piece =
    "begin mrbayes;\n" .
      "set autoclose=yes nowarn=yes;\n";
  my $seed_piece = '';
  if (defined $self->{seed}) {
    my $seed = $self->{seed}; $seed_piece .=  "set seed=$seed;\n";
  } 
  if (defined $self->{swapseed}) {
    my $swapseed = $self->{swapseed}; $seed_piece .= "set swapseed=$swapseed;\n";
  }
  my $middle_piece =  "execute $alignment_nex_filename;\n" .
    "set precision=6;\n" .
      "lset rates=invgamma;\n" .
	"prset aamodelpr=fixed(wag);\n" .
	  "$prset_pinvarpr" . 
	    # "prset pinvarpr=fixed(0.15);\n" .
	    "mcmcp minpartfreq=0.02;\n" . # bipartitions with freq. less than this are not used in the  diagnostics (default is 0.10)
	      "mcmcp allchains=yes;\n" .
		"mcmcp burninfrac=$burnin_frac;\n" .
		  "mcmcp nchains=$n_temperatures;\n" .
		    "mcmcp nruns=$n_runs;\n" .
		      "mcmcp temp=$temperature_gap;\n" .
			"mcmcp samplefreq=$sample_freq;\n" .
			  "mcmcp printfreq=$print_freq;\n" .
			    #  "mcmcp filename=$file_basename;\n" .
			    "mcmcp checkpoint=yes checkfreq=$chunk_size;\n";
  my $end_piece = "sump;\n" . "sumt;\n" . "end;\n";

  $self->{mrbayes_block1} = 
    $begin_piece . $seed_piece .
      $middle_piece . "mcmc ngen=$chunk_size;\n" .
	$end_piece;

  $self->{mrbayes_block2} = 
    $begin_piece . $middle_piece .
      "mcmc append=yes ngen=$chunk_size;\n" .
	$end_piece;

  return $self;
}

sub run{
  my $self = shift;

  my $chunk_size = $self->{chunk_size};
  my $ngen = $chunk_size;
  my $mrbayes_block1 = $self->{mrbayes_block1};

  open my $fh, ">tmp_mrb1.nex";
  print $fh "$mrbayes_block1";
  close $fh;

  my $mb_output_string = `mb tmp_mrb1.nex`;
  $self->{ngens_run} = $ngen;

  my $mc3swap_filename = $self->{file_basename} . ".mc3swap";
  open my $fhmc3, ">$mc3swap_filename";
  print $fhmc3 "$ngen ", $self->extract_swap_info($mb_output_string), "\n";

  open $fh, ">mb1.stdout";
  print $fh "$mb_output_string \n";
  close $fh;

  my ($converged, $conv_string) = $self->test_convergence($self->{file_basename});
  my $converge_count += $converged;
  open my $fhc, ">MB.converge";
  print $fhc "$ngen $converge_count  $conv_string\n";

  foreach (my $i=1; $i>0; $i++) { # infinite loop
    $ngen += $chunk_size;
    my $mrbayes_block2 = $self->{mrbayes_block2};
    $mrbayes_block2 =~ s/ngen=\d+;/ngen=$ngen;/; # subst in the new ngen

    open $fh, ">tmp_mrb2.nex";
    print $fh "$mrbayes_block2";
    close $fh;

    $mb_output_string =  `mb tmp_mrb2.nex`;

    $self->{ngens_run} = $ngen;

    print $fhmc3 "$ngen ", $self->extract_swap_info($mb_output_string), "\n";
    open $fh, ">mb2.stdout";
    print $fh "$mb_output_string \n";
    close $fh;

    ($converged, $conv_string) = $self->test_convergence($self->{file_basename});
    $converge_count += $converged;
    print $fhc "$ngen $converge_count  $conv_string\n";
    last if($converge_count >= $self->{converged_chunks_required});
  }
  close $fhmc3; close $fhc;
  return;
}


sub splits_convergence{
  my $self = shift;
  my $file_basename = shift;	# e.g. fam9877
  my $min_hits = $self->{splits_min_hits}; # ignore splits with fewer hits than this.
  my $max_ok_stddev = $self->{splits_max_ok_stddev}; # convergence is 'ok' for a split if stddev < this.
  my $max_ok_avg_stddev = $self->{splits_max_ok_avg_stddev}; # convergence is 'ok' for a split if stddev < this.
  #print "in splits convergence file basename: ", $file_basename, "\n"; #exit;

  my $filename = $file_basename . ".nex.tstat";
  # print "filename: $filename\n";
  open my $fh, "<$filename";
  my @lines = <$fh>;

  my ($avg_stddev, $count, $bad_count) = (0, 0, 0);
  foreach (@lines) {
    #   print;
    next unless(/^\s*\d/); # skip if first non-whitespace is not numeral.
    my @cols = split(" ", $_);
    my $hits = $cols[1]; my $stddev = $cols[3];
    #  print "$hits, $min_hits, $stddev\n";
    last if($hits < $min_hits);
    $count++;
    $avg_stddev += $stddev;
    if ($stddev > $max_ok_stddev) {
      $bad_count++;
      next;
    }
  }
  $avg_stddev = ($count == 0)? 100 : $avg_stddev/$count;
  my $splits_converged = ($bad_count == 0  and  $avg_stddev < $max_ok_avg_stddev);
  return ($splits_converged, $count, $bad_count, $avg_stddev);
}


sub modelparam_convergence{	# look at numbers in *.pstat file
  # to test convergence
  my $self = shift;
  my $file_basename = shift;
  my $min_ok_ESS = $self->{modelparam_min_ok_ESS};
  my $max_ok_PSRF = $self->{modelparam_max_ok_PSRF};
  my $max_ok_KSD = $self->{modelparam_max_ok_KSD};
  my $string = '';
  my $ngens_skip = int($self->{burnin_frac} * $self->{ngens_run});

  open my $fh, "<$file_basename.nex.pstat";
  my @lines = <$fh>;
  close $fh;
  my $discard = shift @lines;
  my $count_param_lines = 0;
  my $KSD_datacol = 1;
  my $LL_KSD = $self->KSDmax($ngens_skip, 1);
  my $n_bad = ($LL_KSD <= $max_ok_KSD)? 0 : 1;
  my @KSDmaxes = ($LL_KSD);
  #  $n_bad++ if($LL_KSD > $self->{modelparam_max_ok_KSD}); # require LogL just to pass KSD test
  foreach (@lines) {
    my @cols = split(" ", $_);
    my ($avgESS, $PSRF) = @cols[7, 8]; # col 6 is the min ESS among the runs, col 7 is avg.
    next unless($avgESS =~ /^\d*[.]?\d+/ and $PSRF =~ /^\d*[.]?\d+/);
    $KSD_datacol++; # 2,3,4,... the params in pstat file are in cols 2,3,4,... in *.run?.p
    my $KSDmax = $self->KSDmax($ngens_skip, $KSD_datacol);
    push @KSDmaxes, $KSDmax;
    $string .= "$avgESS $PSRF ";
    if ($avgESS < $min_ok_ESS
	or $PSRF > $max_ok_PSRF
	or  $KSDmax > $max_ok_KSD) {
      $n_bad++;
    }
  }
  $string .=  join(" ", map sprintf("%5.3f", $_), @KSDmaxes); #    join(" ", @KSDmaxes);
  return ($n_bad, $string);
}

sub test_convergence{
  my $self = shift;
  my $file_basename = shift;	# e.g. fam9877.nex

  my ($splits_converged, $splits_count, $splits_bad_count, $splits_avg_stddev) =
    $self->splits_convergence($file_basename);
  my ($modelparam_n_bad, $modelparam_string) =
    $self->modelparam_convergence($file_basename);
  my $ngens_skip = int($self->{burnin_frac} * $self->{ngens_run});
 
  my $conv_string = "$splits_count $splits_bad_count $splits_avg_stddev " .
    " $modelparam_string  $modelparam_n_bad  ";
 
  my $converged = ($splits_converged  and  $modelparam_n_bad == 0);
 
  return ($converged? 1 : 0, $conv_string);
}


sub extract_swap_info{
  my $self = shift;
  my $mb_stdout_string = shift;
  my @mb_stdout_lines = split("\n", $mb_stdout_string);
  my $n_lines_to_extract = 0;
  my $extract_next_n = 0;
  my $n_runs = undef;
  my $n_chains = undef;
  my $out_string = '';
  foreach (@mb_stdout_lines) {
    if (/number of chains to (\d+)/) {
      $n_chains = $1;
      $n_lines_to_extract = $n_chains + 4;
      last if(defined $n_runs);
    } elsif (/number of runs to (\d+)/) {
      $n_runs = $1;
      last if(defined $n_chains);
    }

  }
  my $run;
  my %run_string = ();
  foreach (@mb_stdout_lines) {
    if (/Chain swap information for run (\d+)/) {
      $run = $1;
      $extract_next_n = $n_lines_to_extract;
    }
    $out_string .= "$_\n" if($extract_next_n > 0);
    $extract_next_n--;
    if ($extract_next_n == 0) {
      $run_string{$run} = $out_string;
      $out_string = '';
      last if($run == $n_runs);
    }
  }
  $out_string = '';

  foreach (keys %run_string) {
    #print "$_  ", $run_string{$_}, "\n";
    my @lines = split("\n", $run_string{$_});
    splice @lines, 0, 4;
    #  print join("\n", @lines);

    my %ij_swap_pA = ();
    my %ij_swap_tries = ();
    foreach my $i (1..$n_chains) {
      my $l = $lines[$i-1];
      $l =~ s/^\s*\d+\s+[|]\s+//;
      my @xs = split(" ", $l);
      my $n_ntry = $i-1;
      my $n_pA = $n_chains-$i;
      foreach my $j (1..$n_ntry) {
	#	print "swap_tries key: [$i $j]\n";
	$ij_swap_tries{"$i $j"} = shift @xs;
      }
      foreach (1..$n_pA) {
	my $j = $_ + $i;
	#	print "swap_pA key: [$j $i]\n";
	$ij_swap_pA{"$j $i"} = shift @xs;
      }
    }				# loop over chains
    my %ij_swap_accepts = ();
    # my @sijs = sort {$a cmp $b} keys %ij_swap_tries;
    # foreach (@sijs) {

    foreach my $diff (1..$n_chains-1) {
      foreach my $i (1..$n_chains-1) {
	my $j = $i + $diff;

	last if($j > $n_chains);
	my $key = "$j $i";
	#	print "i,j: $i, $j, key: [$key] \n";
	if (exists $ij_swap_pA{$key} and exists $ij_swap_tries{$key}) {
	  $ij_swap_accepts{$key} = $ij_swap_tries{$key} * $ij_swap_pA{$key};
	  $out_string .= int($ij_swap_accepts{$key}+0.5) . " " . $ij_swap_tries{$key} . "  ";
	} else {
	  warn "key $key present in neither ij_swap_tries nor ij_swap_pA.\n";
	}
      }
      $out_string .= ' ';
    }
    $out_string .= ' ';
  }				# loop over runs
  return $out_string;
}


sub KSDmax{	     # This does all pairwise comparisons between runs
  # for one of the params in the *.run?.p files
  # and returns the largest Kolmogorov-Smirnov D
  my $self = shift;

  my $ngen_skip = shift || 0;
  my $datacol = shift;		# data column to use.
  $datacol = 1 if(!defined $datacol);

  my $bigneg = -1e300;

  my $file_basename = $self->{alignment_nex_filename}; # e.g. fam9877
  # store data in hashes
  my @val_count_hrefs = ({}, {}); #
  my @counts = (0, 0);
  my @files = `ls $file_basename.run?.p`;
  my $runs_to_analyze = scalar @files;
  warn "in KSDmax. n_runs: ", $self->{n_runs}, " *.p files found: ", $runs_to_analyze, " should agree, using min of the two.\n" if($self->{n_runs} != $runs_to_analyze);
  $runs_to_analyze = min($runs_to_analyze, $self->{n_runs});

  foreach my $irun (1..$runs_to_analyze) {
    my $i = $irun - 1;
    my $filename = "$file_basename.run" . $irun . ".p"; 
    open my $fh, "<$filename";
    while (<$fh>) {
      my @cols = split(" ", $_);
      # skip non-numerical stuff.
      next unless($cols[0] =~ /^\d+$/); 
      my ($ngens, $x) = @cols[0,$datacol];
      next if($ngens < $ngen_skip);
      $val_count_hrefs[$i]->{$x}++;
      $counts[$i]++;
    }
    close $fh;
  }

  # get cumulative distributions:
  my @val_cumeprob_hrefs = ();
  foreach my $i (0..$runs_to_analyze-1) {
    push @val_cumeprob_hrefs, TlyUtil::cumulative_prob($val_count_hrefs[$i], $counts[$i]);
  }

  my @KSDs = (); # Kolmogorov-Smirnov D for each pairwise comparison between runs
  foreach my $i (0..scalar @files - 2) {
    foreach my $j ($i+1..scalar @files - 1) {
      push @KSDs, TlyUtil::Kolmogorov_Smirnov_D(@val_cumeprob_hrefs[$i,$j]);
    }
  }
  return max(@KSDs);
}

