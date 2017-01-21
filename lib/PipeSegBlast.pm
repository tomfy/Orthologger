package PipeSegBlast;
use strict;
use base 'PipelineSegment';
use List::Util qw (min max sum);

sub new {
   my $class = shift;
   my $control_filename = shift; # 
   my $args  = {};
   my $self  = bless $args, $class;
   return $self;
}
