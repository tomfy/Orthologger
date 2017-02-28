package PipelineSegment;
use strict;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use Cwd qw (getcwd abs_path);

sub new {
   my $class = shift;
   my $pl_obj = shift;
   my $args  = {'attributes' => [], 'hidden' => {}};
   my $self  = bless $args, $class;
   $self->set('pipeline', $pl_obj);
   print STDERR "done constructing PipelineSegment \n";
   return $self;
}

sub init{
   my $self = shift;
   my $defaults = shift;
   print STDERR "In PipelineSegment->init() sets defaults.\n";
   while (my ($k, $v) = each %$defaults) {
      $self->set($k, $v);
   }
}


sub get{
   my $self = shift;
   my $attr = shift;
   return (exists $self->{$attr})?  $self->{$attr} : undef;
}

sub set{
   my $self = shift;
   my $attr = shift;
   my $value = shift;
   if (!exists $self->{$attr}) { # this attribute is new - didn't previously exist.
      if (!exists $self->{attributes}) { # so store it in attributes array
         $self->{attributes} = [$attr];
      } else {
         push @{$self->{attributes}}, $attr;
      }
   }
   $self->{$attr} = $value;
}

sub stringify{
   my $self = shift;
   my $string = '';
   $string .= 'Segment: ' . sprintf("%s \n", (ref $self or $self));
   for my $k (@{$self->get('attributes')}) {
      if (exists $self->get('hidden')->{$k} and $self->get('hidden')->{$k}){ # skip if this attribute is hidden (e.g. gg_string)
      $string .= "$k  [hidden] \n";
      } else {
         my $v = $self->get($k) // 'undef';
         $string .= "  $k  $v \n";
      }
   }
   return $string;
}

1;
