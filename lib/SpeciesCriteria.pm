package SpeciesCriteria;
use strict;
use List::Util qw (min max sum);

# Stores a specification of criteria on a set of species.
# Then, given a set of species, can report whether the criteria
# are satisfied.

sub new {
   my $class = shift;
   my $filename_or_string = shift;
  my $self  = bless {}, $class;
   my $in_string = '';
   if (-f $filename_or_string) {
      open my $fh, "<", "$filename_or_string" or die "couldn't open $filename_or_string for reading.\n";
      while (my $line = <$fh>) {
         $in_string .= $line;
      }
      close $fh;
   }else{
      $in_string = $filename_or_string;
   }
   my %name_speciessetstr = ();
   my $sp_content_criteria_string = '';

   my @lines = split("\n", $in_string);
   for my $line (@lines) {
      # print $line, "\n";
      if ($line =~ /^\s*$/) {
         # empty line - do nothing
      } elsif ($line =~ /[(].*[)]\s*[(].*[)]/) {
 #        print STDERR "content criteria line: $line\n";
         $sp_content_criteria_string .= $line . "\n";
      } else {
  #       print STDERR "species set line: $line \n";
         $line =~ /^\s*(\S+)\s+(.*)$/;
         my ($name, $speciessetstr) = ($1, $2);
         $name_speciessetstr{$name} = $speciessetstr;
      }
   }

   my %name_spset = ();
   while (my($name, $spset) = each(%name_speciessetstr)) {
      $name_spset{$name} = {};
      $spset =~ s/\s//g;        #remove whitespace
      my @species = split(',', $spset);
      for my $sp (@species) {
         $name_spset{$name}->{$sp} = 1;
      }
   }
   $self->{name_spset} = \%name_spset;
   $self->{spcontent_criteria} = $sp_content_criteria_string;

   return $self;
}

sub get_speciesset{ #
   my $self = shift;
   my $name = shift;
   if(exists $self->{name_spset}->{$name}){
      return $self->{name_spset}->{$name}
   }else{
      return undef;
   }
}

sub check_criteria{
   my $self = shift;
   my @species_present = @_;
   my %name_spcount = ();
   my $name_spset = $self->{name_spset};
   for my $sp (@species_present) {
      while ( my ($name, $spset) = each %$name_spset) {
         if (exists $spset->{$sp}) {
            $name_spcount{$name}++;
         }
      }
   }

   my @criterion_lines = split("\n", $self->{spcontent_criteria});
   my $OK = 0;
   for my $criterion_line (@criterion_lines) {
   #   print STDERR $criterion_line, "\n";
      my @anded_criteria = split(";", $criterion_line);
      my $line_OK = 1;
      for my $andedcrit (@anded_criteria) {
         if (! $self->check_criterion($andedcrit, \%name_spcount)) {
            $line_OK = 0;
            last;
         }
      }
 #     print STDERR $criterion_line, "  OK?  $line_OK \n";
      if ($line_OK) {
         $OK = 1;
         last;
      }
   }
   return $OK;
}

sub check_criterion{
   my $self = shift;
   my $criterion = shift;       # e.g. (AMd+Ambtr)(10,1000)
   my $spset_count = shift;
   my $count = 0;
   my ($cats, $minmax) = split('\)\(', $criterion);
   $cats =~ s/[(]//;
   $minmax =~ s/[)]//;
   my ($min,$max) = split(',', $minmax);
   my @categories = split('\+', $cats);
   for (@categories) {
      if (exists $spset_count->{$_}) {
         $count += $spset_count->{$_};
      }
   }
   my $OK = 1;
   $OK = 0 if($count < $min);
   $OK = 0 if($count > $max);
 #  print "  category/ies: $cats ; minmax: $minmax ; count: $count ; OK: $OK \n";

   return $OK;
}


1;
