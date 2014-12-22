package Phyml;
use strict;
use List::Util qw (min max sum);

# a object to run phyml (phylogeny inference package)
# phyml takes a phylip-format alignment as input
# which allows only 10 characters for the sequence name!
# much too short. This package allows supplying a fasta file
# or string to the constructor, the sequence names are then replaced
# with seq_0, seq_1, etc. and later, after the tree is inferred 
# these are replaced with the original (correct) sequence names.


my $default_arg_href = {
			'dna_or_protein' => 'protein',
			'interleaved_or_sequential' => 'sequential',
			'fasta_string' => undef,
			'fasta_file' => undef,
			'phylip_string' => undef,
			'phylip_file' => undef,
			'n_bootstrap' => 0,
			'optimize_param' => 'none', # l, lt, ltr
			'initial_tree_newick_file' => undef,
			'alpha' => undef,
			'n_rate_categories' => undef,
			'p_invariant' => undef,
			'subst_model' => undef
		       };


sub new {
	my $class = shift;
	my $arg_href = shift;

	my $self = bless {}, $class;

# initialize parameters to defaults.
	foreach (keys %$default_arg_href) {
		$self->{$_} = $default_arg_href->{$_};
	}
# reset any parameters specified in argument hash ref.
	foreach (keys %$arg_href) {
		$self->{$_} = $arg_href->{$_};
	}

	if(defined $self->{fasta_string}){
	  $self->fasta_string_to_phylip($self->{fasta_string});
	  $self->write_phylip_file();
	}elsif(defined $self->{fasta_file}){
		$self->fasta_file_to_phylip(); 
	}else{
		die "Must specify alignment fasta string or fasta file in Phyml constructor argument hash ref.\n";
	}

	$self->construct_phyml_command_line();


} # end of constructor

sub construct_phyml_command_line{
	my $self = shift;

	my $phyml_command_line = 'phyml ';
	if($self->{dna_or_protein} eq 'protein'){ $phyml_command_line .= ' -d aa '; }
	if($self->{interleaved_or_sequential} eq 'sequential'){ $phyml_command_line .= ' -q '; }

	if(defined $self->{phylip_file}){ $phyml_command_line .= ' -i ' . $self->{phylip_file}; }
	if(defined $self->{n_bootstrap}){ $phyml_command_line .= ' -b ' . $self->{n_bootstrap}; }
	if(defined $self->{optimize_param}){ $phyml_command_line .= ' -o ' . $self->{optimize_param}; }
	if(defined $self->{alpha}){ $phyml_command_line .= ' -a ' . $self->{alpha}; }
	if(defined $self->{n_rate_classes}){ $phyml_command_line .= ' -c ' . $self->{n_rate_classes}; }
	if(defined $self->{p_invariant}){ $phyml_command_line .= ' -v ' . $self->{p_invariant}; }
	if(defined $self->{subst_model}){ $phyml_command_line .= ' -m ' . $self->{subst_model}; }

	if (defined $self->{initial_tree_newick}) {
	  my $init_tree_newick = $self->{initial_tree_newick};
	  my $init_newick_encoded = $self->encode_newick($init_tree_newick);
#	  print "encoded init newick: $init_newick_encoded\n";

	  my $tmp_filename = $$ . "_ft.newick_tmp";
	  open my $fhinittree, ">","$tmp_filename";
	  print $fhinittree "$init_newick_encoded;\n";
	  $phyml_command_line .= " -u $tmp_filename --quiet ";
	}

# 	if(defined $self->{initial_tree_newick_file}){
# 	  my $init_tree_newick_file = $self->{initial_tree_newick_file};
# 	  open my $fhinittree, "<", "$init_tree_newick_file";
# 	  my $init_newick = join("", <$fhinittree>);
# #print "init newick $init_newick\n";
	  
#   $init_newick = $self->encode_newick($init_newick);
# #print "encoded init newick: $init_newick\n";
# 	  close $fhinittree; 


# open $fhinittree, ">","init_newick_tmp";
# 	  print $fhinittree "$init_newick\n";
# 	  $phyml_command_line .= ' -u init_newick_tmp --quiet '; 
# 	}

#	print STDERR "phyml command line: $phyml_command_line \n";
	$self->{phyml_command_line} = $phyml_command_line;
	return $self;
}

sub fasta_file_to_phylip{
# get the input alignment into proper format:
# takes a fasta alignment file name as argument,
# gets phylip format string and writes to file (and stores string and filename in obj.)
	my $self = shift;
	my $new_fasta_file = shift || undef;

	my $align_fasta_file;
	if(defined $new_fasta_file){
		$align_fasta_file = $new_fasta_file;
	}elsif(defined $self->{fasta_file}){
		$align_fasta_file = $self->{fasta_file};
	}else{
		die "In Phyml::fasta_file_to_phylip. No fasta file specified.\n";
	}

# set up "encoded" ids of form seq_n (n=0,1,2...) and 1-1 mapping between 
# actual sequence ids and these
	my %id_encodedid = ();
	my %encodedid_id = ();
	my %encodedid_seq = ();

	my $seq_number = 0;
	my $seq_length;

	open my $fh, "<", "$align_fasta_file";
	while(<$fh>){

		if(/^>(\S+)/){
			my $id = $1;
			my $encoded_id = 'seq_' . $seq_number;
			$id_encodedid{$id} = $encoded_id;
			$encodedid_id{$encoded_id} = $id;
			my $seq = <$fh>;
			$seq =~ s/^\s*(\S.*\S)\s*$/$1/; # remove initial, final whitespace.

				if($seq_number > 0){
					warn "Sequences in alignment have different lengths.\n" if(length $seq ne $seq_length);
				}
			$seq_length = length $seq;

			$encodedid_seq{$encoded_id} = $seq;
			$seq_number++;
		}else{
			print STDERR "In Phyml::fasta_file_to_phylip. No initial > ; line skipped: $_";
		}
	}
	close $fh;

	$self->{id_encodedid} = \%id_encodedid;
	$self->{encodedid_id} = \%encodedid_id;
	$self->{encodedid_seq} = \%encodedid_seq;

# get phylip formatted string and print to file
	my $phylip_string = "$seq_number $seq_length \n";
	my $line_length = 80;
	foreach (keys %encodedid_seq){
		my $id_and_seq = substr($_ . '           ', 0, 10) . $encodedid_seq{$_};
		my $pos = 0;
		while($pos < length $id_and_seq){
			$phylip_string .= substr($id_and_seq, $pos, $line_length) . "\n";
			$pos += $line_length;
		}
	}
	my $phylip_file = $align_fasta_file;
	$phylip_file =~ s/fasta$//;
	$phylip_file =~ s/[.]$//;
	$phylip_file .= '.sqphlp';
	open $fh, ">", "$phylip_file";
	print $fh $phylip_string;
	close $fh;

	$self->{phylip_file} = $phylip_file;
	$self->{phylip_string} = $phylip_string;
#return $phylip_file;

}

sub fasta_string_to_phylip{
# get the input alignment into proper format:
# takes a fasta alignment string as argument,
# gets phylip format string and stores it in obj.
	my $self = shift;
	my $fasta_string = shift || undef;

# set up "encoded" ids of form seq_n (n=0,1,2...) and 1-1 mapping between 
# actual sequence ids and these
	my %id_encodedid = ();
	my %encodedid_id = ();
	my %encodedid_seq = ();

	my $seq_number = 0;
	my $seq_length;

my @fasta_lines = split("\n", $fasta_string);
#	while(<$fh>){
	while (@fasta_lines) {
	  $_ = shift @fasta_lines;
	  if (/^>(\S+)/) {
	    my $id = $1;
	    my $encoded_id = 'seq_' . $seq_number;
	    $id_encodedid{$id} = $encoded_id;
	    $encodedid_id{$encoded_id} = $id;
	    my $seq = shift @fasta_lines;
	      $seq =~ s/^\s*(\S.*\S)\s*$/$1/; # remove initial, final whitespace.
	    if ($seq_number > 0) {
	      warn "Sequences in alignment have different lengths.\n" if(length $seq ne $seq_length);
	    }
	    $seq_length = length $seq;

	    $encodedid_seq{$encoded_id} = $seq;
	    $seq_number++;
	  } else {
	    print STDERR "In Phyml::fasta_string_to_phylip. No initial > ; line skipped: $_";
	  }
	}
#	close $fh;

	$self->{id_encodedid} = \%id_encodedid;
	$self->{encodedid_id} = \%encodedid_id;
	$self->{encodedid_seq} = \%encodedid_seq;

# get phylip formatted string and print to file
	my $phylip_string = "$seq_number $seq_length \n";
	my $line_length = 80;
	foreach (keys %encodedid_seq){
	  my $sequence $encodedid_seq{$_};
	  $sequence =~ s/[J]/-/g; # Phyml can't handle amino-acid 'J' (leucine/isoleucine); replace with - (gap)
		my $id_and_seq = substr($_ . '           ', 0, 10) . $sequence; # key (id) should be <= 10 char

		my $pos = 0;
		while($pos < length $id_and_seq){
			$phylip_string .= substr($id_and_seq, $pos, $line_length) . "\n";
			$pos += $line_length;
		}
	}
	$self->{phylip_string} = $phylip_string;
	return $phylip_string;
      }

sub write_phylip_file{
  my $self = shift;
  my $phylip_file;
  if (defined $self->{phylip_file}) {
    $phylip_file = $self->{phylip_file};
 
  } elsif (defined $self->{fasta_file}) { # if phylip filename not already defined, construct from fasta file name.
    $phylip_file = $self->{fasta_file};
    $phylip_file =~ s/fasta$//;
    $phylip_file =~ s/[.]$//;
    $phylip_file .= '.sqphlp';
    $self->{phylip_file} = $phylip_file;
  } else {
    my $PID = $$;
    $phylip_file = $PID . "_align.phylip";
  }
  $self->{phylip_file} = $phylip_file;
#  print STDERR "phylip filename: $phylip_file \n";
  open my $fh, ">", "$phylip_file" or die "Couldnt open $phylip_file for reading. \n";
  if ($self->{phylip_string}) {
    print $fh $self->{phylip_string};
  } else {
    die "Phylip string: [", $self->{phylip_string}, "]??? \n";
  }
  close $fh;
}

sub run{
  my $self = shift;
  my $phyml_command_line = $self->{phyml_command_line};


  #print "Phyml command line: \n";
  #print "$phyml_command_line \n";
  #	system "$phyml_command_line";
  my $phyml_stdout = `$phyml_command_line`;
  $self->{phyml_stdout} = $phyml_stdout;
  my $lnl = ($phyml_stdout =~ /Log likelihood of the current tree:\s+(\S+)/)? $1 : '---';
  my $cpu_time = ($phyml_stdout =~ /Time used (\d+)h(\d+)m(\d+)s/)? 3600*$1 + 60*$2 + $3 : '---';
  $lnl =~ s/[.]\s*$//;		# remove final . if present.
  #print "lnL: $lnl \n";
  $self->{log_likelihood} = $lnl; # store likelihood
$self->{cpu_time} = $cpu_time;
  my $phyml_tree_outfile = $self->{phylip_file};
  #	$phyml_tree_outfile =~ s/seqphylip$//; 
  my $phyml_stats_outfile = $phyml_tree_outfile;
  $phyml_tree_outfile .= "_phyml_tree.txt";
  $phyml_stats_outfile .= "_phyml_stats.txt";

  my $newick = `cat $phyml_tree_outfile`;
  if (0) {
    my @sorted_encodedids = sort { length $b <=> length $a } keys %{$self->{encodedid_id}};

    foreach my $encoded_id (@sorted_encodedids) {
      #       print $enc_id, "\n";
      my $id = $self->{encodedid_id}->{$encoded_id};
      $newick =~ s/$encoded_id/$id/g;
    }
  } else {
    $newick = $self->decode_newick($newick);
  }

  #print $newick, "\n";

  $newick =~ s/,/,\n/g;
  #print $newick, "\n";
  $self->{newick_out} = $newick;
  # print STDERR "in Phyml::run. returning.\n";
}

sub encode_newick{
  my $self = shift;
  my $newick = shift;
	my @sorted_ids = sort { length $b <=> length $a } keys %{$self->{id_encodedid}};

	foreach my $id (@sorted_ids){
		my $encoded_id = $self->{id_encodedid}->{$id};
	#	print "id, encid: $id  $encoded_id\n";
	#	$newick =~ s/$id/$encoded_id/g;
		# can't use regex here because of special characters in
		# sequence names (e.g. | .
	$newick = literal_substitution($newick, $id, $encoded_id);
	}
  return $newick;
}

sub decode_newick{
  my $self = shift;
  my $newick = shift;
	my @sorted_encodedids = sort { length $b <=> length $a } keys %{$self->{encodedid_id}};

	foreach my $encoded_id (@sorted_encodedids){
#       print $enc_id, "\n";
		my $id = $self->{encodedid_id}->{$encoded_id};
		$newick =~ s/$encoded_id/$id/g;
	}
  return $newick;
}

sub literal_substitution{
  my $string = shift;
  my $oldpart = shift;
my $replacementpart = shift;

my $oldlength = length $oldpart;
  foreach (0..((length $string) - $oldlength)){
    if(substr($string, $_, $oldlength) eq $oldpart){
      substr($string, $_, $oldlength, $replacementpart);
    }
  }
return $string;
}

sub get_phyml_command_line{
  my $self = shift;
  if(defined $self->{phyml_command_line}){
    return $self->{phyml_command_line};
  }else{
    return 'undefined';
  }
}

1; 
