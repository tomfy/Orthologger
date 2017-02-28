package PipeSegFamilies;
use strict;
use base 'PipelineSegment';
use List::Util qw (min max sum);
use Cwd qw (getcwd abs_path);

sub new {
   my $class = shift;
   my $pl_obj = shift;          #
   #   my $input_filenames = shift; # whitespace separated .m8 filenames
   my $args  = {                #  defaults 
                output_dir => 'families',
                size_limit => '500',
                max_e_value => 1e-6,
                multiplicity_knee => 1000,
                log10_eval_penalty => 0,
                #  fastas_output_dir => 'family_fastas',
                input_filename_pattern => '../*.m8', # by default (if no predecessor segment output file defined), look in ../
               };
   my $self  = bless $args, $class;
   $self->set('pipeline', $pl_obj); # store the pipeline obj to which this segment belongs.
   #   $self->set('input_filenames', $input_filenames);
   return $self;
}

sub run{
   my $self = shift;

   my $pl_obj = $self->get('pipeline');
   my $gg_filename =  $pl_obj->get('gg_filename');
   my $fh_progress = $pl_obj->get('fh_progress');

   my $family_max_e_value = $self->get('max_e_value');
   my $family_size_limit = $self->get('size_limit');
   my $family_multiplicity_knee = $self->get('multiplicity_knee');
   my $family_log10_eval_penalty = $self->get('log10_eval_penalty');

   print STDERR "top of PipeSegFamilies->run; ", getcwd(), " ", abs_path(getcwd()), "\n";
   my $fam_abcs_dir;
   if (defined $self->get('predecessor')) {
      print STDERR "predecessor:  ", $self->get('predecessor'), "\n";
      $fam_abcs_dir = $self->get('output_dir');
      mkdir $fam_abcs_dir unless(-d $fam_abcs_dir); # create dir if doesn't exist
      chdir($fam_abcs_dir);
   } else { # should be in dir (created by hand) which is subdir of blast output dir. should have *.m8 files
      $fam_abcs_dir = getcwd();
   }

   # *************** make family (abc) files *********************
   my $m8_filename_str = `ls ../*.m8`;
   die "ls ../*.m8 finds no files.\n" if($m8_filename_str =~ /No such file or directory/);

   my @m8_filenames = split(" ", $m8_filename_str);
   my @abc_part_filenames = ();
   for my $m8_filename (@m8_filenames) { # $self->get('input_filenames'))) {
      print STDERR "m8 filename: $m8_filename \n";
      my ($v, $dir, $fname) = File::Spec->splitpath($m8_filename);
      $fname =~ s/m8$/abc/;

      my $fam_abcs_filename = $fname; # $fam_abcs_dir . '/' . $fname;
      my $m8tofams_cl = "m8toabc_fams.pl -input $m8_filename -max_eval $family_max_e_value -gg $gg_filename -output $fam_abcs_filename ";
      $m8tofams_cl .= " -fam_size_limit $family_size_limit ";
      $m8tofams_cl .= "-multiplicity_knee $family_multiplicity_knee -log10penalty $family_log10_eval_penalty ";

      print STDERR "m8tofams_cl: $m8tofams_cl \n";
      # exit;
      my $abc_filename = `$m8tofams_cl`;
      # -fam_size_limit $family_size_limit `;
      $abc_filename =~ s/\s+$//; 
      print STDERR "abc_filename: [$abc_filename]\n";
      push @abc_part_filenames, $abc_filename;
   }
   #my $abc_part_filenames = `split_abc.pl $abc_filename $n_pieces`;
   #my @abc_parts = split(" ", $abc_part_filenames);
   #exit;

   print "[" ,join("][", @abc_part_filenames), "]\n";
   print $fh_progress "abc family filenames: \n", join("\n", @abc_part_filenames), "\n";

   # ************** done making family (abc) files *******************
   die "exiting after making family (abc)files.\n";
   # ************* make family (fastas)

   my $family_taxon_requirement = $pl_obj->get('family_taxon_requirement');
   my $added_groups_string = $pl_obj->get('added_groups_string');

   my $all_species_fasta_filename = $pl_obj->get('all_species_fasta_filename');
   die "fasta input file undefined or file does not exist.\n" unless(defined $all_species_fasta_filename and -f $all_species_fasta_filename);
   my $fam_fastas_dir = $pl_obj->get('families_output_dir');
   mkdir $fam_fastas_dir unless(-d $fam_fastas_dir);
   my @fastas_filenames = ();
   for my $the_abc_part (@abc_part_filenames) {
      my ($v, $dir, $fname) = File::Spec->splitpath($the_abc_part);
      $fname =~ s/abc$/fastas/;
      my $output_fastas_filename = $fam_fastas_dir . '/' . $fname;
      #   print STDERR "ZZZZ: ", `which seq+matches2fasta.pl`, "\n";
      #  my $addgrpstring = "'" . 'amposdicots:Solanum_lycopersicum,Theobroma_cacao,Medicago_truncatula,Manihot_esculenta;monocots:Oryza_sativa,Phoenix_dactylifera' . "'";
      my $seqm2f_cl = "seq+matches2fasta.pl -gg $gg_filename -abc_file $the_abc_part -fasta_infile $all_species_fasta_filename -output_filename $output_fastas_filename";
      $seqm2f_cl .= " -taxon_requirement $family_taxon_requirement " if(defined $family_taxon_requirement);
      $seqm2f_cl .= " -added_groups $added_groups_string " if(defined $added_groups_string);
      print STDERR "seq+... cl:  $seqm2f_cl \n";
      system("$seqm2f_cl");
      # "seq+matches2fasta.pl -gg $gg_filename -abc_file $the_abc_part -fasta_infile $all_species_fasta_filename -output_filename $output_fastas_filename");
      push @fastas_filenames, $output_fastas_filename
   }

   my $fastas_filenames_str = join("\n", @fastas_filenames);
   $self->set('fastas_filenames', $fastas_filenames_str);
   print STDERR "Family fastas file ready to align: \n", $fastas_filenames_str, "\n";

   print $fh_progress "Family_fastas_files: " , $fastas_filenames_str, "\n";
}                               # end of make_families section





1;
