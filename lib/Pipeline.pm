package Pipeline;
use strict;
use List::Util qw (min max sum);

sub new {
   my $class = shift;
   my $control_filename = shift; # 
   my $defaults = shift;
   my $args  = {};
   my $self  = bless $args, $class;
   $self->initialize($control_filename, $defaults);

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

sub get_hash_val{
   my $self = shift;
   my $href = shift;
   my $key = shift;
   if ( (exists $self->{$href}) and 
        (ref $self->{$href} eq 'HASH') and 
        (exists $href->{$key}) ) {
      return $href->{$key};
   } else {
      return undef;
   }
}

sub set_hash_val{
   my $self = shift;
   my $href = shift;
   my $key = shift;
   my $value = shift;

   if ( (exists $self->{$href}) and
        (ref $self->{$href} eq 'HASH') ) {
      $self->{$href}->{$key} = $value;
   } else {
      $self->{$href} = {$key => $value};
   }
}

1;


