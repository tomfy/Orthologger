package PipelineSegment;
use strict;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use Cwd qw (getcwd abs_path);
use File::Spec qw(splitpath);

use TomfyMisc qw(_stringify);

use parent 'OrderedHash';

sub new {
   my $class = shift;
   print STDERR "Top of PipelineSegment constructor.\n";
   my $self = $class->SUPER::new(); # OrderedHash obj.
   $self->merge(
                segment_name => undef,
                pipeline => undef, # put these in with undef values here so order of keys is set.
                state => undef,
                predecessor => undef,
                successor => undef,
                output_dir => undef,
                # hidden => OrderedHash->new(), # pipeline => 1),
               );
   $self->hide('pipeline');
   return $self;
}

sub init{
   my $self = shift;
   $self->merge(@_);
   $self->set(state => 'initialized'); # state can be: initialized, complete
}

sub stringify{
   my $self = shift;
   my $indent_increment = shift || '  ';
   my $spacer1 = shift || '  '; # between key and value
   my $spacer2 = shift || "\n"; # between key/value pairs of hash, or between elements of array
   my $string = '';
   my $indent = '  ';
   $string .= "segment " . (ref $self) . "  " . $self->get('segment_name') . "\n";
   $string .= "{\n";
   for my $a ($self->keys()) {
      my $v = $self->get($a) // 'undef';
      if ( $self->_hidden($a)) {
         #  $string .= "$indent$a$spacer1" . "hidden$spacer2"; # . (ref $v or $v) . ",\n";
      } else {
         $string .= "$indent$a$spacer1" . _stringify($v, $indent . $indent_increment, $spacer1, $spacer2) . "$spacer2";
      }
   }
   $string .= "}\n";
   return $string;
}

sub get_predecessor{
   my $self = shift;
   my $predecessor_seg_name = $self->get('predecessor');
   if (defined $predecessor_seg_name) {
      return $self->get('segments')->get($predecessor_seg_name);
   } else {
      return undef;
   }
   return $predecessor_seg_name;
}

sub set_successor_state{ # set the state of the successor segment to $new_state 
   # if there is no successor segment, do nothing.
   my $self = shift;
   my $new_state = shift;
   my $successor_name = $self->get('successor');
   if (defined $successor_name) {
      my $successor_obj = $self->get('pipeline')->get('segments')->get($successor_name);
      $successor_obj->set(state => $new_state);
   }
}

sub awaken{
   my $self = shift;
   # do nothing - override in derived classes for segments to actually do something.
   # idea is to recover the object from its output, instead of rerunning everything
   # e.g. for PipeSegGG, get the seqid_species hash from the output file.
   return;
}

1;
