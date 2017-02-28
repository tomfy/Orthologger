package Pipeline;
use strict;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use 5.012;                      # so can use each with arrays

use parent 'Hash::Ordered';

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
   print STDERR "Top of Pipeline::initialize() \n";

   my $added_groups_string = '';
   my ($seg, $seg_class, $seg_obj) = ('pipeline', undef, undef);
   my $segments = Hash::Ordered->new();

   open my $fh_ctrl, "<", "$control_filename" or die "Could not open $control_filename for reading; exiting.";
   while (my $line = <$fh_ctrl>) {
      next if($line =~ /^\s*#/); # skip all comment lines
      next if($line =~ /^\s*$/); # skip lines that are all whitespace
      $line =~ s/#.*$//;         # delete comment at end of line.
      if ($line =~ /^\s*__END__/) {
         last;
      } elsif ($line =~ /^\s*segment\s+(\S+)\s+(\S+)/) {
         ($seg, $seg_class) = ($1, $2);
         print STDERR "About to construct  $seg  $seg_class \n";
         my $evalstring = 'print STDERR "AAA\n"; $seg_obj = ' . "$seg_class" . '->new($self); ';
         print STDERR "About to 'eval' expression: $evalstring \n";
         my $evalval = eval $evalstring;
         print STDERR "evalval: [$evalval] \n";
         print STDERR (defined $evalval)? $evalval : 'undef', "\n";
         print STDERR "Constructed: ", ref $seg_obj, "\n";
         $segments->set($seg => $seg_obj);
         print STDERR "ZZZ: [", $seg_obj->stringify(), "]\n";
      } else {
         if ($seg eq 'pipeline') {
            my @strs = split(" ", $line);
            if (scalar @strs == 2) {
               #     $self->{$strs[0]} = $strs[1];
               $self->set($strs[0] => $strs[1]);
            } elsif (scalar @strs == 3) {
               print STDERR "ABC: ", join("; ", @strs), "\n";
               if (!defined $self->get($strs[0])) {
                  #   $self->{$strs[0]} = Hash::Ordered->new($strs[1] => $strs[2]);
                  $self->set($strs[0] => Hash::Ordered->new($strs[1] => $strs[2]));
               } else {
                  #    $self->{$strs[0]}->set($strs[1] => $strs[2]);
                  $self->get($strs[0])->set($strs[1] => $strs[2]);
               }
            } else {
               warn "Unexpected line in control_file: $line";
            }
         } else {
            if ($line =~ /^\s*(\S+)\s+(\S+)/) {
               print STDERR "setting PipelineSegment key/obj pair to $1 $2 \n";
               $seg_obj->set($1 => $2);
            }
         }
      }
   }
   my $include = 0;
   my $included_segments = Hash::Ordered->new();
   for ($segments->keys()) {
      $include = 1 if($_ eq $self->get('start'));
      $included_segments->set($_ => $segments->get($_)) if($include);
      $include = 0 if($_ eq $self->get('last'));
   }
   $segments = $included_segments;
   my @segnames = $segments->keys();
   for my $i (keys @segnames) {
      my $seg_obj = $segments->get($segnames[$i]);
      my $predecessor_obj = ($i>0)? $segments->get($segnames[$i-1]) : undef;
      print STDERR "$i  $segnames[$i] ", 
        (blessed $seg_obj)? blessed $seg_obj : 'undef', "  ", 
          (blessed $predecessor_obj)? blessed $predecessor_obj : 'undef', "\n"; # exit;
      $seg_obj->set('predecessor' => $predecessor_obj);
   }
   $self->set('segments' => $segments);
   print STDERR $self->stringify(), "\n";
   return;
}

sub stringify{
   my $self = shift;
   my $string = "Pipeline: " . sprintf("%s \n", (ref $self or $self));


   for my $k ($self->keys()) {
      my $v = $self->get($k);
      $string .= "$k" if($k ne 'segments');
      if (blessed($v) and $v->isa('Hash::Ordered')) {
         $string .= "\n";
         for my $ik ($v->keys()) {
            my $iv = $v->get($ik);
            if (blessed($iv) and $iv->isa('PipelineSegment')) {
               #  $string .= $iv->stringify() . "\n\n";
            } else {
               $string .= "  $ik  $iv\n";
            }
         }
      } elsif (ref $v eq 'HASH') {
         $string .= "\n";
         while (my($ik, $iv) = each %$v) {
            $string .= "  $ik  " . (ref $iv || $iv) . "\n";
         }
      } elsif (ref $v eq 'ARRAY') {
         $string .= "\n  ";
         $string .= join("\n  ", @$v) . "\n";
      } else {
         $string .= "  $v\n";
      }
   }

   for my $seg_obj ($self->get('segments')->values()) {
      $string .= "\n" . $seg_obj->stringify() . "\n";
   }
   $string .= "\n ** end ** \n";
   return $string;
}


1;


