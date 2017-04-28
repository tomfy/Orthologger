package TomfyMisc;
use Exporter qw 'import';
@EXPORT_OK = qw 'run_quicktree run_fasttree run_phyml store_gg_info timestring date_time file_format_is_abc_iie split_fasta_file split_fasta_string fix_fasta_files short_species_name read_block clean _stringify _destringify';
use strict;
use Scalar::Util qw(blessed);

sub timestring{
   my $s_time = shift;          # time in seconds (e.g. from time() )
   my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
   my @ltime = localtime($s_time);
   my $t_string = $ltime[2] . ":" . $ltime[1] . ":" . $ltime[0];
   my $d_string = $ltime[5]+1900 . "-" . $months[$ltime[4]] . "-" . $ltime[3];
   return "$t_string, $d_string.";
}

sub date_time{
   my @ltime = split(" ", localtime);
   my $date = join('', @ltime[4,1,2]); # year month day e.g. 2017Jan21
   my $time = $ltime[3];               # hour:min:sec

   return ($date, $time);
}

sub store_gg_info {	 #xx    # read in gene-genome association file
   # store in hash (keys: seqids, values: species). e.g. ('ATxxxx' => Arabidopsis_thaliana');
   my $gg_filename   = shift;
   my %seqid_species = ();
my %species_count = ();
   if ( defined $gg_filename ) {
      if ( -f $gg_filename ) {
         open my $fh_gg, "<", "$gg_filename" or die "Couldn't open $gg_filename for reading.\n";
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
                  $species_count{$species}++;
               }
            }
         }
         close $fh_gg; # done storing gg_file info in hash %seqid_species
      } else { # 
         die "$gg_filename: no such file.\n";
      }
   } else {
      die "gg filename undefined. \n";
   }
   return (\%seqid_species, \%species_count);
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

sub split_fasta_file{
   my $fasta_filename = shift;
   my $n_parts = shift;

   #   my ($v, $dir, $fname) = File::Spec->splitpath($fasta_filename);
   #   print STDERR "[$v], [$dir], [$fname] \n";
   my $wcout = `grep '^>' $fasta_filename | wc`;
   $wcout =~ /^\s*(\S+)/;
   my $n_seqs = $1;
   my $n_sequences_in_each_part = int($n_seqs/$n_parts) + 1;

   open my $fh, "<", "$fasta_filename" or warn "couldn't open $fasta_filename for reading.\n";

   my $target_sequence_count = $n_sequences_in_each_part;
   print "total number of sequences: $n_seqs ; number of parts: $n_parts approx. number of sequences per part: $target_sequence_count \n";

   my $i_part = 1;
   my $output_filename_stem = $fasta_filename;
   $output_filename_stem =~ s/[.]fasta$//;
   my $output_filename = $output_filename_stem . "_part$i_part.fasta";
   my @output_fasta_part_filenames = ($output_filename);
   open my $fh_out, ">", $output_filename;
   my $sequence_count = 0;
   while (<$fh>) {
      if (/^>/) {
         s/^(>\s*\S+).*/$1/;    # just keep > and id.
         $sequence_count++;
         #       print "$sequence_count  $_\n";
         if ($sequence_count > $target_sequence_count) {
            close $fh_out;
            $i_part++;
            $output_filename = $output_filename_stem . "_part$i_part.fasta";
            $target_sequence_count += $n_sequences_in_each_part;
            push @output_fasta_part_filenames, $output_filename;
            open $fh_out, ">", $output_filename;
         }
      }
      print $fh_out  $_;
   }
   close $fh_out;
   close $fh;
   return \@output_fasta_part_filenames;
}

sub split_fasta_string{
   my $fasta_string = shift;
   my $n_parts = shift;

   my @lines = split("\n", $fasta_string);
   my $n_lines = scalar @lines;
   my $lines_so_far_this_part = 0;
   my $lines_per_part = 1 + int ($n_lines/$n_parts);
   my @parts = ();
   my $part_string = '';
   for my $a_line (@lines) {
      #   print "n parts: $n_parts  n_lines: $n_lines lines per part: $lines_per_part  $lines_so_far_this_part \n";
      if ($a_line =~ /^>/) {
         if ($lines_so_far_this_part >= $lines_per_part) {
            print STDERR "lines in this part: $lines_so_far_this_part.\n";

            #   print "\n\n", $part_string, "\n\n";
            push @parts, $part_string;
            $part_string = '';
            $lines_so_far_this_part = 0;
         }
      }
      $part_string .= $a_line . "\n"; #print $part_string;
      $lines_so_far_this_part++;
   }
   # print "\n\n", $part_string, "\n\n";
   push @parts, $part_string;
   print STDERR "lines in this part: $lines_so_far_this_part. RETURNING FROM split_fasta_string.\n";
   return \@parts;
}

sub fix_fasta_files{            #
   my $min_seqlength = shift;
   my @fasta_files = @_;
   my $fasta_output_string = '';
   # my $min_seqlength = $self->get('min_sequence_length');
   my ($count, $short_count) = (0, 0);

   for my $fasta_file (@fasta_files) {
      open my $fh, "<", glob($fasta_file) or die "Couldn't open $fasta_file for reading.\n";
      my $idline;
      while ($idline = <$fh>) {
         if ($idline =~ /^>\s*(\S+)/) {
            $idline = ">$1\n";
            last;
         }
      }
      my $sequence = '';
      while (my $the_line = <$fh>) {
         if ($the_line =~ /^>\s*(\S+)/) {
            my $new_idline = ">$1\n"; # remove any stuff after id. 
            $count++;
            if (length $sequence >= $min_seqlength) { # keep only sequences which are long enough.
               $fasta_output_string .= $idline . $sequence . "\n";
            } else {
               $short_count++;
            }
            if ($count % 100000 == 0) {
               print STDERR "$count  $short_count \n";
            }
            $sequence = ''; # done with one sequence; reset to start next sequence.
            $idline =  $new_idline;
         } else {                     # sequence line
            $the_line =~ s/[*\s]+$//; # remove final * and whitespaces
            $sequence .= $the_line;
         }

      }
      if (length $sequence >= $min_seqlength) {
         $fasta_output_string .= $idline . $sequence . "\n";
      } else {
         $short_count++;
      }
      close $fh;
      print STDERR "all sequences: $count  short sequences: $short_count \n";

   }                            # end of loop over files
   return $fasta_output_string;
}

sub short_species_name{ # return abbreviated form of name: First letter of genus and 2 first letters of species
   my $genus_species = shift;
   if ($genus_species =~ /^(\S)\S+_(\S{2})/) {
      return $1 . $2;
   } else {
      return $genus_species;
   }
}

sub read_block{                # read a block defined by { and }
   # other {} blocks can be nested inside.
   # returns a string containing the block, starting with {, ending with }
   my $fh = shift;

   my $string = '';
   my $curly_count = 0;
   $_ = <$fh>;
   if (s/^\s*{/{/) {
      $curly_count += tr/{// - tr/}//;
      $string .= $_;
   } else {
      die "In read_block, first line: $_, is not start of block.\n";
   }
   while (<$fh>) {
      $string .= $_;
      $curly_count += tr/{// - tr/}//;
      last if($curly_count == 0);
   }
   print "read_block return value [$string]\n";
   return $string;
}

sub clean{                      # remove comments and blank lines.
   my $string = shift;
   my @lines = split("\n", $string);
   $string = '';
   for (@lines) {
      s/#.*$//;                 #delete from first # to end of line
      $string .= "$_\n" if(/\S/); # include lines with at least one non-whitespace character.
   }
   return $string;
}

# sub decomment{ # for each line, remove from first # to end of line.
#    my $string = shift;
#    my @lines = split("\n", $string);
#    for (@lines) {
#       s/#.*$//;                 #delete from first # to end of line
#    }
#    return join("\n", @lines);
# }

# sub dewhitespace{ # remove lines which contain only whitespace
#    my $string = shift;
#    my @lines = split("\n", $string);
#    $string = '';
#    for (@lines) {
#       $string .= "$_\n" if(/\S/); # include lines with at least one non-whitespace character.
#    }
#    return $string;
# }

sub _destringify{
   my $string = shift;
   my $obj = shift || undef;
   my $result;
#   print STDERR "TOP of _destringify. string: [|$string|} \n"; #
   if ($string =~ s/^\s*{//) {  # hash (Hash::Ordered)
      my $curly_count = 1;
      $result = (defined $obj)? $obj : Hash::Ordered->new();
      die "In {} block. obj should be a Hash::Ordered, but is ", ref $result, ".\n" if(! $result->isa('Hash::Ordered'));
      while ($curly_count > 0) {
#         print STDERR "in _destringify; string: $string \n";
         if ($string =~ s/^\s*[}]//) {
            $string =~ s/^\s*,//;
            $curly_count--;
         } elsif ($string =~ s/^\s*([^\s=>]+)(\s*=>\s*|\s+)([^\s=>]+)//) {
            my $key = $1;
            my $value = $3;
            $value =~ s/,\s*$//;  # remove final comma if present
            $string =~ s/^\s*,//; #remove initial comma if present
            if ($value =~ /^\s*[{[]/) {
               $string = $value . $string;
               ($value, $string) = _destringify($string);
            }
            $result->set($key, $value);
         } else {
            die "problem in _destringify {}: $string \n";
         }
      }  # end while loop
   } elsif ($string =~ s/^\s*[[]//) { # array ref.
      my $square_count = 1;
      die "In [] block. Object should be undefined or array ref, but is ", ref $obj, ".\n" if(defined $obj and (ref $obj ne 'ARRAY'));
      $result = [];
      while ($square_count > 0) {
         if ($string =~ s/^\s*(\S+)//) {
            my $value = $1;
            $value =~ s/,\s*$//;  # remove final comma if present
            $string =~ s/^\s*,//; # remove initial comma if present
            if ($value =~ /^\s*[]]/) {
               $square_count--;
            } else {
               if ($value =~ /^\s*[{[]/) {
                  $string = $value . $string;
                  ($value, $string) = _destringify($string);
               }
               push @$result, $value;
            }
         } else {
            die "problem in _destringify []: $string \n";
         }
      }
   }
#   print STDERR "In _destringify, [$string] \n"; #exit;
   return ($result, $string);
}

sub _stringify{
   my $x = shift;      #  Hash::Ordered, hashref, arrayref, or scalar;
   my $indent = shift;
   my $spacer1 = shift || ' ';  # between key and value
   my $spacer2 = shift || "\n"; # between key/value pairs of hash, or between elements of array
   my $str = '';
   if ((blessed $x) and $x->isa('Hash::Ordered')) {
  #    print STDERR "hash::ordered branch.\n";
      $str .= "{\n";
      for my $k ($x->keys()) {
         my $v = $x->get($k);
         $str .= $indent . $k . $spacer1 . _stringify($v, $indent . '  ', $spacer1, $spacer2) . "$spacer2"; # "$spacer2\n";
      }
      $str =~ s/$spacer2$/\n$indent}/;
   } elsif (ref $x eq 'HASH') {
  #    print STDERR "hashref branch.\n";
      $str .= "{\n";
      for my $k (keys %$x) {
         my $v = $x->{$k};
         my $xxx = $indent . $k . $spacer1 . _stringify($v, $indent . '  ', $spacer1, $spacer2) . "$spacer2"; # \n";
         $str.= $xxx; # $indent . $k . $spacer1 . _stringify($v, $indent . '  ', $spacer1, $spacer2) . "$spacer2\n";
     #    print STDERR "v: $v, xxx: $xxx";
      }
      $str =~ s/$spacer2$/\n$indent}/;
   } elsif (ref $x eq 'ARRAY') {
   #   print STDERR "arrayref branch.\n";
      $str .= "[\n";
      for my $v (@$x) {
         $str .= $indent . _stringify($v, $indent . '  ', $spacer1, $spacer2) . "$spacer2";
      }
      $str =~ s/$spacer2$/\n$indent]/;
   } elsif ((ref $x) eq '') {   # not a ref
   #   print STDERR "not a ref branch. x: ", $x // 'undef', "\n";
      $str .= (defined $x)? sprintf("%s", $x) : 'undef' . "\n";
   } else {
   #   print STDERR "other ref branch\n";
      $str .= ref $x . "\n";
   }
   return $str;
}


1;
