package OrderedHash;
use strict;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use 5.012;                      # so can use each with arrays

use TomfyMisc qw (_destringify);

use parent 'Hash::Ordered';

sub new {
   my $class = shift;
   my $self = $class->SUPER::new();  # Hash::Ordered obj.

   # $self->merge(
   #              segment_name => undef,
   #              pipeline => undef, # put these in with undef values here so order of keys is set.
   #              state => undef,
   #              predecessor => undef,
   #              successor => undef,
   #              output_dir => undef,
   #              hidden => {pipeline => 1},
   #             );
   print STDERR "done constructing OrderedHash \n";
   return $self;
}

sub hide{                 # hide all the items in the argument list @_
   my $self = shift;
my $hidden = $self->get('hidden');
   if(defined $hidden){
      die "in OrderedHash::hide; hidden is not a Hash::Ordered.\n" if(! $hidden->isa('Hash::Ordered'));
   }else{
      $self->set(hidden => Hash::Ordered->new());
   }
   for (@_) {
      $self->get('hidden')->set($_ => 1);
   }
}

sub _hidden{
   my $self = shift;
   my $a = shift;
   my $hidden = $self->get('hidden');
   if(defined $hidden){
      return $hidden->get($a) // 0;
   }else{
      return 0;
   }
}

1;
