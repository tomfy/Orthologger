#!/usr/bin/perl
use strict;
use List::Util qw ( min max sum);
use Math::GammaFunction qw ( log_gamma );
use Math::Random qw ( random_multinomial );
use Math::GSL::SF qw ( gsl_sf_lngamma );

# my $sample_size = shift || 100;

#my $nstring = shift || '3,5,7,10,2,5,7';
my $pstring = shift || '1,2,4,8,16,32';
my $Ndata = 100;
my $multiplier = 3;
print "#ps: $pstring\n";
my $ref_sample_size = 1000;

#my @ns = split(",", $nstring);
my @ps = split( ",", $pstring );
my $sump = sum(@ps);
for (@ps) {
  $_ /= $sump;
  print "#ps: ", join( "; ", @ps ), "\n";
}

for (1..1000) {
  my %NNNN_count = ();
  my @data_sample = random_multinomial( $Ndata, @ps );
  my $Nincrement = int($Ndata/20);
  my $old_count = -100;
  for (my $NNNN = 12*$Ndata; $NNNN > 0; $NNNN -= $Nincrement ) {
 
    my @data_x = @data_sample;
    for (@data_x) {
      $_ *= $NNNN/$Ndata; # scale so sum of counts is $NNNN
    }
    my $log_prob_data_x = log_multinomial_prob( \@data_x, \@ps);

    my $count = 0;
#    my @log_probs = ();
    for ( 1 .. $ref_sample_size ) {
      my @drawn_ns = random_multinomial( $NNNN, @ps );
      my $log_prob = log_multinomial_prob( \@drawn_ns, \@ps );
#      push @log_probs, $log_prob;
      $count++ if($log_prob_data_x > $log_prob);
    }
    if($count >= $ref_sample_size/2   and  $old_count < $ref_sample_size/2){
      my $answer = $NNNN + 0.5*$Nincrement;
print "eff sample size: $answer \n";
last;
    }
    $old_count = $count;
    # @log_probs = sort { $a <=> $b } @log_probs;
    # my $mid = int( $ref_sample_size / 2 );
    # my $median_log_prob = ($log_probs[$mid] + $log_probs[$mid+1])/2;
    # $NNNN_count{$NNNN} = $count;

    # #unshift @log_probs, -1e300;
    # #push @log_probs, 1e300;

    # my $min_log_prob = $log_probs[1] - 0.001;
    # my $max_log_prob = $log_probs[ scalar @log_probs - 2 ] + 0.001;
  }
  # my @NNNNs = sort {$b <=> $a} keys %NNNN_count;
  # for (@NNNNs) {
  #   if ( $NNNN_count{$_} > $ref_sample_size/2) {
  #     print "eff sample size: $_ \n";
  #     last;
  #   }
  # }
}

sub log_multinomial_prob{ # when elements of @ns are integers (>=0) should give same
  # as log_multinomial_prob, but also works for non-integers in @ns.
  my @ns = @{ (shift) };
  my @ps = @{ (shift) };
  my $N = sum(@ns);
  my $sump = sum(@ps);		# normalize probabilities
  for (@ps) {
    $_ /= $sump;
  }

  my $result = log_gamma($N+1);
  while (my ($i, $mi) = each @ns) {
    my $pi = $ps[$i];
    $result += $mi*log($pi) 
      - log_gamma($mi+1);	# function from Math::GammaFunction
    #   - gsl_sf_lngamma($mi+1); # same thing using gsl
  }
  return $result;
}


sub binary_locate {
  my $value  = shift;
  my @sarray = @{ (shift) };
  my $lo     = shift || 0;
  my $hi     = shift || scalar @sarray - 1;
  die "value $value not in range [$lo, $hi]\n"
    if ( $value < $sarray[$lo] or $value > $sarray[$hi] );
  my $count= 0;
  while ( ( $hi - $lo ) > 1 ) {
    my $i = int( ( $hi + $lo ) / 2 );
    if ( $value < $sarray[$i] ) {
      $hi = $i;
    } else {
      $lo = $i;
    }
    $count++;
  }
  return ( $lo, $hi, $count );
}


