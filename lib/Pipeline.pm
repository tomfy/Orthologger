package Pipeline;
use strict;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use Cwd qw (getcwd abs_path);
use 5.012;                      # so can use each with arrays

use TomfyMisc qw(date_time short_species_name clean read_block _stringify _destringify);

use parent 'OrderedHash';       # 'Hash::Ordered';

sub new {
   my $class = shift;
   my $self = $class->SUPER::new(); # $params->as_list());
   print STDERR "in Pipeline constructor. class arg, class of obj constructed : $class   $self \n";
   return $self;
}

sub initialize{
   my $self = shift;
   my $control_filename = shift;
   my $params = shift;

   $self->merge($params->as_list());
   my $xxx = '';
   my ($pipe_class, $pipe_name) = (undef, undef);
   my ($seg_class, $seg_name, $seg_obj) = (undef, undef, undef);
   my $segments = Hash::Ordered->new();
   open my $fh_ctrl, "<", "$control_filename" or die "Could not open $control_filename for reading; exiting.";
   while (my $line = <$fh_ctrl>) {
      next if($line =~ /^\s*#/); # skip all comment lines
      next if($line =~ /^\s*$/); # skip lines that are all whitespace
      $line =~ s/#.*$//;         # delete comment at end of line.
      last if($line =~ /^\s*__END__/);
      if ($line =~ /^\s*pipeline\s+(\S+)\s+(\S+)?/) {
         ($pipe_class, $pipe_name) = ($1, $2); # Not used at present March 2017
         $self->set('pipeline_name' => $pipe_name);
         my $pipeline_block = clean(read_block($fh_ctrl));
         ($self, $xxx) = _destringify($pipeline_block, $self);
         print STDERR "Done destringifying pipeline block.\n";
      } elsif ($line =~ /^\s*segment\s+(\S+)\s+(\S+)?/) {
         ($seg_class, $seg_name) = ($1, $2);
         my $evalstring = 'print STDERR "AAA\n"; $seg_obj = ' . "$seg_class" . '->new($seg_name, $self); ';
         print "About to 'eval' expression: $evalstring \n";
         my $evalval = eval $evalstring;
         print "evaleval: ", (defined $evalval)? "[$evalval]" : 'undef', "\n";
         die "evaleval is undefined.\n" if(!defined $evalval);
         print "Constructed: ", ref $seg_obj, "\n";
         my $segment_block = clean(read_block($fh_ctrl));
         ($seg_obj, $xxx) = _destringify($segment_block, $seg_obj);
         print "ref seg_obj: ", ref $seg_obj, "\n";
         $segments->set($seg_name => $seg_obj);
      } elsif (!defined $pipe_name ) {
         die "Pipeline name is undefined. Bye.\n";
      }
   }

   # set up predecessor and successor segments
   my @segnames = $segments->keys();
   for my $i (keys @segnames) {
      my $seg_obj = $segments->get($segnames[$i]);
      print STDERR "i, segname: $i  ", $segnames[$i], "\n";
      $seg_obj->set(predecessor => ($i>0)? $segnames[$i-1]: undef);
      $seg_obj->set(successor => ($i< (scalar @segnames-1))? $segnames[$i+1] : undef);
      print STDERR "Predecessor:  ", $seg_obj->get('predecessor') // 'undef', "  successor: ", $seg_obj->get('successor') // 'undef', "\n";
   }
   # # set up 'next' state
   # for my $i (keys @segnames) {
   #    my $seg_obj = $segments->get($segnames[$i]);
   #    my $pred_seg_name = $seg_obj->get('predecessor');
   #    if (defined $pred_seg_name) {
   #       my $pred_seg_obj = $segments->get('$pred_seg_name');
   #       if ($pred_seg_obj->get('state') eq 'complete' and $seg_obj->get('state') ne 'complete') {
   #          $seg_obj->set(state => 'next');
   #          print STDERR "A  ", $pred_seg_name // 'undef  ', $segnames[$i], "  ", $seg_obj->get('state'), "\n";
   #          last;
   #       }
   #    } else {                  # no predecessor 
   #       if ($seg_obj->get('state') ne 'complete') {
   #          $seg_obj->set(state => 'next');
   #          print STDERR "B  ", $pred_seg_name // "undef  ", $segnames[$i], "  ", $seg_obj->get('state'), "\n";
   #          last;
   #       }
   #    }
   #    exit;
   # }

   # find the first segment which is not 'complete', and set its state to 'next'.
   # for my $i (keys @segnames) { 
   #    my $seg_obj = $segments->get($segnames[$i]);
   #    my $seg_state = $seg_obj->get('state');
   #    if (!defined $seg_state  or  $seg_state ne 'complete') {
   #       $seg_obj->set(state => 'next');
   #       last;
   #    }
   # }

   $self->set('segments' => $segments);
   print STDERR "bottom of Pipeline->initialize. \n";
   return;
}

sub run{
   my $self = shift;
   my $segments = $self->get('segments');
   my $iter = $segments->iterator();
   while (my ($segment_name, $segment_obj) = $iter->()) {

      print STDERR "segment name and state: ", $segment_obj->get('segment_name'), "  ", $segment_obj->get('state'), "\n";

      if ($segment_obj->get('state') eq 'complete') {
         print STDERR "about to call awaken method of ", $segment_name, " ", $segment_obj, "\n";
         $segment_obj->awaken();
         print STDERR $segment_obj->get('segment_name'), " awakened.\n"; #exit;
      } else {
         print STDERR "about to call run method of ", $segment_name, " ", $segment_obj, "\n";
         $segment_obj->run();
         print STDERR "after run method of ", $segment_name, " ", $segment_obj, " wd is: ", abs_path('./'),  "\n";

         my $state_filename = $self->get('state_filename');
         open my $fh, ">", "$state_filename" or die "couldn't open $state_filename for writing.\n";
         print $fh $self->stringify();
      }
   }
}

sub stringify{
   my $self = shift;
   my $indent_increment = shift || '  ';
   my $spacer1 = shift || '  '; # between key and value
   my $spacer2 = shift || "\n"; # between key/value pairs of hash, or between elements of array
   my $string = '';
   my $indent = '  ';
   $string .= "pipeline " . (ref $self) . "  " . $self->get('pipeline_name') . "\n";
   $string .= "{\n";
   for my $a ($self->keys()) {
      next if($a eq 'segments');
      my $v = $self->get($a) // 'undef';
      if ( $self->_hidden($a)) { # suppress output
         # $string .= "$indent$a$spacer1" . "hidden$spacer2"; # . (ref $v or $v) . ",\n";
      } else {
         #  print STDERR "GGGGG a: $a v: $v \n";
         $string .= "$indent$a$spacer1" . _stringify($v, $indent . $indent_increment, $spacer1, $spacer2) . "$spacer2";
      }
   }
   $string .= "}\n";            # end of pipeline block.
   for my $seg_obj ($self->get('segments')->values()) {
      $string .= $seg_obj->stringify() . "\n";
   }
   #   $string .= "}\n";
   #  $string .= "** end ** \n";

   return $string;
}

sub reset_start_output_dir{ # resets the first not 'complete' segment's output dir
   # to $new_output_dir 
   my $self = shift;
   my $new_output_dir = shift;
   my $segments = $self->get('segments');
   for my $segment_name ($segments->keys()) {
      my $segment_obj = $segments->get($segment_name);
      if ($segment_obj->get('state') ne 'complete') {
         $segment_obj->set(output_dir => abs_path($new_output_dir));
         print STDERR "Segment $segment_name has just had its output dir reset to: ", abs_path($new_output_dir), "\n"; #exit;
         last;
      }
   }
}


1;




