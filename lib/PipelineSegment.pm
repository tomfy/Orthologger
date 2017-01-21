package PipelineSegment;
use strict;
use List::Util qw (min max sum);

sub new {
   my $class = shift;
#   my $control_filename = shift; # 
   my $args  = {};
   my $self  = bless $args, $class;
return $self;
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
   $self->{$attr} = $value;
}

1;
