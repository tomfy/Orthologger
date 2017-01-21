package PhylogenomicPipeline;
use strict;
use List::Util qw (min max sum);
use base 'Pipeline';

sub new {
   my $class = shift;
   my $control_filename = shift; # 
   my $defaults = shift;
   my $args  = {};
   my $self  = bless $args, $class;
   $self->initialize($control_filename, $defaults);

   return $self;
}

sub initialize{
   my $self = shift;
   my $control_filename = shift;
   my $defaults = shift;
   while (my($k,$v) = each %$defaults) {
      $self->{$k} = $v;
   }
   open my $fh_ctrl, "<", "$control_filename" or die "Could not open $control_filename for reading; exiting.";
   my $added_groups_string = '';
   $self->{taxon_file} = {};
   while (<$fh_ctrl>) {
      next if(/^\s*#/);         # skip all comment lines
      next if(/^\s*$/);         # skip lines that are all whitespace
      s/#.*$//;                 # delete comment at end of line.
      my @strs = split(" ", $_);
      if (scalar @strs == 2) {
         $self->{$strs[0]} = $strs[1];
      } elsif (scalar @strs == 3) {
         if (!exists $self->{$strs[0]}) {
            $self->{$strs[0]} = {$strs[1] => $strs[2]};
         } else {
            $self->{$strs[0]}->{$strs[1]} = $strs[2];
         }
      } else {
         warn "Unexpected line in control_file: $_";
      }
   }
   $self->set('n_species', scalar keys %{$self->get('taxon_inputpath')});
   #   print STDERR 'Added groups String: ' . "$added_groups_string \n";
   $added_groups_string =~ s/;\s*$//; # remove final ;
   $added_groups_string = "'" . $added_groups_string . "'";
   $self->{added_groups_string} = $added_groups_string;
   #print STDERR 'Added groups String: ' . "$added_groups_string \n";
   #$self->{taxon_file} = \%taxon_file;
   #  $self->{param_name_val} = \%param_name_val;
   return;
}


sub stringify{
   my $self = shift;
   my $string = "Stringified pipeline object: \n";
   #my @sks = sort {if(ref $self->{$
   my %hrefs = ();
   my %sclrs = ();
   my %arefs = ();
   while (my($k, $v) = each %$self) {
      if (ref $v eq 'HASH') {
         $hrefs{$k} = $v;
      } elsif (ref $v eq 'ARRAY'){ 
         $arefs{$k} = $v;
} else {
         $sclrs{$k} = $v;
      }
   }
   my @skeys = sort keys %arefs;
    for my $akey (@skeys) {
      $string .= "$akey \n" . join("\n  ", @{$arefs{$akey}}) . "\n";
   }
@skeys = sort keys %hrefs;
   for my $akey (@skeys) {
      $string .= "$akey\n";
      while (my($ik, $iv) = each %{$hrefs{$akey}}) {
         $string .= "  $ik  $iv\n";
      }
   }
   @skeys = sort keys %sclrs;
   for my $akey (@skeys) {
      $string .= "$akey " . $sclrs{$akey} . "\n";
   }
   $string .= "\n ** end ** \n";
   return $string;
}

1;


