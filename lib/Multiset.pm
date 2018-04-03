package Multiset;
use Scalar::Util qw (blessed);
use List::Util qw (min max sum);
use 5.012;                      # so can use each with arrays

sub new {
   my $class = shift;
   my $items = shift; 
   my $args = {};
   my $self  = bless $args, $class;
   $self->{'items'} = {};
   $self->populate($items) if($items);
   print STDERR "done constructing Multiset \n";
   return $self;
}

sub populate{
   my $self = shift;
   my $items = shift;           # array ref

   for my $i (@$items) {
      #  print "adding $i to multiset.\n";
      $self->{items}->{$i}++;
   }
   $self->print();
}

sub sum{
   my $self = shift;
   my $other = shift;
   return $self->plus($other);
}

sub plus{
   my $self = shift;
   my $other = shift;
   my $sum = {};

   while ( my ($i,$n) = each %{$self->{items}} ) {
      #   print "x: $i $n\n";
      $sum->{$i} += $n;
   }
   while ( my ($i,$n) = each %{$other->{items}} ) {
      # print "y: $i $n\n";
      $sum->{$i} += $n;
   }
   my $sum_set = Multiset->new();
   #  print "sum items:  ", join(", ", keys %$sum), "\n";
   $sum_set->{items} = $sum;
   return $sum_set;
}

sub diff{
   my $self = shift;
   my $other = shift;
   return $self->minus($other);
}

sub minus{
   my $self = shift;
   my $other = shift;
   my $sum = {};

   while ( my ($i,$n) = each %{$self->{items}} ) {
      #   print "x: $i $n\n";
      $sum->{$i} += $n;
   }
   while ( my ($i,$n) = each %{$other->{items}} ) {
      # print "y: $i $n\n";
      $sum->{$i} -= $n;
   }
   my $sum_set = Multiset->new();
   #  print "sum items:  ", join(", ", keys %$sum), "\n";
   $sum_set->{items} = $sum;
   return $sum_set;
}

sub print{
   my $self = shift;
   my %items_count_hash = %{$self->{items}};
   # print scalar keys %items_count_hash, "\n";
   while (my($i, $n) = each %{$self->{items}}) {
      print "$i  $n \n";
   }
}

  1;
