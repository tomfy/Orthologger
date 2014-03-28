#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );
use Math::Symbolic;
use Algorithm::CurveFit;
use Data::Dumper;

my $n_runs = shift;
my $n_temps = shift;
my $dT = shift || 1;
my @pas = ();
for(1..$n_temps){
	push @pas, [(0) x ($n_runs+1)];
}

my @deltats = ();
my @paccepts = ();


my @lines = <>;
my $last_line = $lines[-1];
my $end_burn_in_line = $lines[ int( (scalar @lines )/10)];
my @ll = split(" ", $last_line);
my @ebil = split(" ", $end_burn_in_line);
my $lgen = shift @ll;
my $ebigen = shift @ebil;
my @pbi_cols = ();
for (0..scalar @ll -1){
	my $diff = $ll[$_] - $ebil[$_];
	push @pbi_cols, $diff;
}

while(@pbi_cols){

	for my $i_run (1..$n_runs){
		for my $gap (1..$n_temps-1){
			my $n_pairs = $n_temps - $gap;
			for my $t_lo (0..$n_pairs-1){
				my $t_hi = $t_lo + $gap;
				my $n_accepted = shift @pbi_cols;
				my $n_out_of = shift @pbi_cols; 
				my $p_accept = $n_accepted / $n_out_of;
				if($t_lo == 0){		
	# print "$lgen $i_run $t_lo $t_hi $p_accept \n";
					$pas[$t_hi]->[$i_run] = $p_accept;
		push @deltats, $t_hi*$dT;
	push @paccepts, $p_accept;
		}
			}
		}
	}
}
for my $t_hi (1..$n_temps-1){
printf("%3i %7.4f  ", $t_hi, $t_hi * $dT);
my $sum = 0;
for my $i_run (1..$n_runs){
my $the_pa = $pas[$t_hi]->[$i_run];
$sum += $the_pa;
printf("%7.5f ", $pas[$t_hi]->[$i_run]);
} 
printf("%7.5f \n", $sum/$n_runs);
}


my ($x1, $y1) = ($deltats[0], $paccepts[0]);
my ($x2, $y2) = ($deltats[-1], $paccepts[-1]);
my $s = (log($y2) - log($y1))/($x2 - $x1);
my $c_guess = $y1/exp($x1*$s);
my $a_guess = $s;
print "#  initial a, c: $a_guess, $c_guess \n";
my $formula = 'c * exp(a*x)';
my $variable = 'x';
my @parameters = (['a', $a_guess, '0.00001'], ['c', $c_guess, '0.00001']);
my $max_iter = 100;
my $square_residual = Algorithm::CurveFit->curve_fit(
	formula => $formula,
	params => \@parameters,
	variable => $variable,
	xdata => \@deltats,
	ydata => \@paccepts,
	maximum_iterations => $max_iter
);
my $a_best = $parameters[0]->[1];
my $c_best = $parameters[1]->[1];
print "#  best a, c: $a_best, $c_best \n";
print "#  rms residual: ", ($square_residual/($n_runs*($n_temps-1)))**0.5, "\n";

my ($dt_best, $pa_best, $max_mstj) = (-1, -1, -1);
$dt_best = -2/$a_best;
$pa_best = $c_best * exp($a_best * $dt_best);
$max_mstj = $pa_best * $dt_best**2;
printf("#  Opt. Tgap, Pa, msTjump: %8.5f %8.5f %8.5f \n", $dt_best, $pa_best, $max_mstj);
