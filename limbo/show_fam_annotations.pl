#!/usr/bin/perl -w
use strict;

#usage e.g.: 
#  show_fam_annotations.pl ORTHOMCL5_ < all_orthomcl.out 
# optional argument is a pattern to match, to select desired cluster,
# otherwise do all of them in all_orthomcl.out file

my $pattern = shift || undef;
my $max_n_descript = shift || 3;

my %taxon_taxonpattern = (     'arabidopsis' => '^AT',
			       'castorbean' => '\d{5}[.]m\d{6}',
			       'medicago' => '^IMGA[|](Medtr|AC|CU)',
			       'rice' => 'LOC_Os',
			       'tomato' => 'Solyc'
			 );
my %taxon_annotfile = (
		       'arabidopsis' => '~/Genome_data/Arabidopsis/Tair10/Tair10_all_pep.fasta',
		       'castorbean' => '~/Genome_data/Castorbean/TIGR_castorWGS_release_0.1.aa.fasta',
		       'medicago' => '~/Genome_data/Medicago/Mt3.5_GenesProteinSeq_20100825.fasta',
		       'rice' => '~/Genome_data/Rice/TIGR6.1/TIGR6.1_pep_nonTE.fasta',
		       'tomato' => '~/Genome_data/Tomato/ITAG2.3/ITAG2.3_proteins.fasta'
		      );

my @taxa = ('arabidopsis', 'castorbean', 'medicago', 'rice', 'tomato');

while (my $line = <>) {
  my %taxon_descripts = (
			 'arabidopsis' => {},
			 'castorbean' => {},  
			 'medicago' => {}, 
			 'rice' => {}, 
			 'tomato' => {}
			);


  #print "line, pattern: ", substr($line, 0, 50), " $pattern \n";;
  #print "match? [", ($line =~ /$pattern/), "]\n";
  next if((defined $pattern) and ! ($line =~ /$pattern/));
  $line =~ s/(ORTHO[^:]*:)\s*//;
  my $cluster_id = $1;
  my @ids = split(" ", $line);
  print "Cluster: $cluster_id \n";
  foreach my $id (@ids) {
    foreach my $taxon (keys %taxon_taxonpattern) {
      my $taxon_pattern = $taxon_taxonpattern{$taxon};
      if ($id =~ /$taxon_pattern/) {
	my $grep_file = $taxon_annotfile{$taxon};
	$id =~ /^([^(]+)[(]\S*[)]/;
	my $grep_out = `grep \'$1\' $grep_file`;
	if ($grep_out =~ /\S/) {
#	  printf("  %10s   ", $taxon);
	  $grep_out =~ /^(\S+)\s*(.*)/;
	  my $description = $2;
	  if ($taxon eq 'tomato') {
	    if ( $grep_out =~ /(["].*["])/ ) {
	      $description = $1;
	    }
	  }
	  $description =~ s/\s+/ /g; 
	  $taxon_descripts{$taxon}->{$description}++;
#	  print "$description \n"; # "$grep_out";
	}
	next;
      }
    }
  }
  for my $txn (@taxa) {
    my $descript_count = $taxon_descripts{$txn};
    my @sort_descripts = sort {length $b <=> length $a } keys %{$descript_count};
    my $old_l = -1;
    my $n_descript_shown = 0;
    foreach (@sort_descripts) {
      if(length $_ ne $old_l){
	print $txn, "  ", $_, "\n" if(length $_ ne $old_l);
      $old_l = length $_;
	$n_descript_shown++;
	last if ($n_descript_shown == $max_n_descript);
      }
    }
  }
  print "\n\n";
}

