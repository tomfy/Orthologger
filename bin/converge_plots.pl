#!/usr/bin/perl -w
use strict;
use lib '/usr/share/perl/5.14.2/';
use Graphics::GnuplotIF qw(GnuplotIF);
my $base_filename = shift || 'famX';
my $do_the_plots  = shift || 1;
my $enhanced      = 1;
my $persist       = 0;

open my $fh0, "<first_chunk.stdout";
my $n_runs;
my $n_temps;
while (<$fh0>) {
  $n_temps = $1 if(/Setting number of chains to (\d+)/);
  if (/Setting number of runs to (\d+)/) {
    $n_runs = $1;
    last;
  }
}
print "n_runs, n_temps: $n_runs, $n_temps \n";

if ( defined $do_the_plots ) {
  my $converge_filename = "$base_filename" . ".converge";
  print "converge filename: $converge_filename \n";
  open (my $fh, "<", $converge_filename) || die "couldn't open $converge_filename for reading.\n";

  my @lines = <$fh>;
  shift
    @lines; # throw away first line - incomplete information for first chunk.

  my $n_cols = scalar split( " ", $lines[0] );

  #print "n cols: $n_cols\n";
  my @col_data = ();
  for ( 1 .. $n_cols ) {
    push @col_data, [];
  }

  while ( my ( $iline, $line ) = each @lines ) {
    my @cols = split( " ", $line );
    while ( my ( $col, $v ) = each @cols ) {

      #print "line: $iline;  col: $col, v: $v.  ref(): [", ref($col_data[$col]), "]\n";
      push @{ $col_data[$col] }, $v;
    }
  }

  #print "\n\n", join(", ", @{$col_data[4]}), "\n";

  my $gens               = $col_data[0];

#  my $topo_L1s           = $col_data[1];
#  my $topo_max_diff = $col_data[2];

  my $splits_avg_stddevs = $col_data[1];
  my $splits_L1s = $col_data[2];
  my $splits_max_diff = $col_data[3];

  my @TL_invESSs     = map( 1 / $_, @{ $col_data[4] } );
  my @alpha_invESSs  = map( 1 / $_, @{ $col_data[5] } );
  my @pinvar_invESSs = map( 1 / $_, @{ $col_data[6] } );

  my $LnL_KSDs    = $col_data[7];
  my $TL_KSDs     = $col_data[8];
  my $alpha_KSDs  = $col_data[9];
  my $pinvar_KSDs = $col_data[10];

  my $LnL_L1s    = $col_data[11];
  my $TL_L1s     = $col_data[12];
  my $alpha_L1s  = $col_data[13];
  my $pinvar_L1s = $col_data[14];

  my $xtics_labels_string =
    " '' 10^2 0, '10^3' 1000 0, '10^4' 10000 0, "
      . " '10^5' 100000 0, '10^6' 1000000 0, '10^7' 10000000 0, "
	. " '' 200 1, '' 300 1, '' 400 1, '' 500 1, '' 600 1, '' 700 1, '' 800 1, '' 900 1, "
	  . " '' 2000 1, '' 3000 1, '' 4000 1, '' 5000 1, '' 6000 1, '' 7000 1, '' 8000 1, '' 9000 1, "
	    . " '' 20000 1, '' 30000 1, '' 40000 1, '' 50000 1, '' 60000 1, '' 70000 1, '' 80000 1, '' 90000 1, "
	      . " '' 200000 1, '' 300000 1, '' 400000 1, '' 500000 1, '' 600000 1, '' 700000 1, '' 800000 1, '' 900000 1, "
		. " '' 2000000 1, '' 3000000 1, '' 4000000 1, '' 5000000 1, '' 6000000 1, '' 7000000 1, '' 8000000 1, '' 9000000 1, "
		  . " '' 20000000 1, '' 30000000 1, '' 40000000 1, '' 50000000 1, '' 60000000 1, '' 70000000 1, '' 80000000 1, '' 90000000 1 ";
  my $no_xtics_labels_string =
    " '' 100 0, '' 1000 0, '' 10000 0, '' 100000 0, '' 1000000 0, '' 10000000 0, "
      . " '' 200 1, '' 300 1, '' 400 1, '' 500 1, '' 600 1, '' 700 1, '' 800 1, '' 900 1, "
	. " '' 2000 1, '' 3000 1, '' 4000 1, '' 5000 1, '' 6000 1, '' 7000 1, '' 8000 1, '' 9000 1, "
	  . " '' 20000 1, '' 30000 1, '' 40000 1, '' 50000 1, '' 60000 1, '' 70000 1, '' 80000 1, '' 90000 1, "
	    . " '' 200000 1, '' 300000 1, '' 400000 1, '' 500000 1, '' 600000 1, '' 700000 1, '' 800000 1, '' 900000 1, "
	      . " '' 2000000 1, '' 3000000 1, '' 4000000 1, '' 5000000 1, '' 6000000 1, '' 7000000 1, '' 8000000 1, '' 9000000 1, "
		. " '' 20000000 1, '' 30000000 1, '' 40000000 1, '' 50000000 1, '' 60000000 1, '' 70000000 1, '' 80000000 1, '' 90000000 1 ";
  my $no_ytics_labels_string =
    " '' 0.01 0, '' 0.1 0, '' 1 0, "
      . " '' 0.003 1, '' 0.004 1, '' 0.005 1, '' 0.006 1, '' 0.007 1, '' 0.008 1, '' 0.009 1, "
	. " '' 0.02 1, '' 0.03 1, '' 0.04 1, '' 0.05 1, '' 0.06 1, '' 0.07 1, '' 0.08 1, '' 0.09 1, "
	  . " '' 0.2 1, '' 0.3 1, '' 0.4 1, '' 0.5 1, '' 0.6 1, '' 0.7 1, '' 0.8 1, '' 0.9 1 ";


  my $scale = 0.5;
  my $plot1 =    # various convergence diagnostic quantities as function of generations
    Graphics::GnuplotIF->new( persist => $persist, style => 'lines lw 2');
  #   $plot1->gnuplot_cmd(' set terminal pdf ');
  #   $plot1->gnuplot_cmd(' set output "converge.pdf" ');
  $plot1->gnuplot_cmd(' set terminal x11 enhanced size 800,600') if ($enhanced);
  $plot1->gnuplot_cmd(' set log ');

  $plot1->gnuplot_cmd(' set multiplot ');
  $plot1->gnuplot_cmd(" set size $scale,$scale ");

# upper two plots
  # splits, topology L1
  #$plot1->gnuplot_set_xrange(1000,100000);
# $plot1->gnuplot_set_title("topological convergence");
  $plot1->gnuplot_cmd(" set xrange [1000:*] ");
  $plot1->gnuplot_cmd(" set yrange [0.002:1] ");
  $plot1->gnuplot_cmd(" set xtics ( $no_xtics_labels_string ) ");
  $plot1->gnuplot_cmd(" set origin 0.0,$scale ");
  $plot1->gnuplot_set_plot_titles( 'splits avg stddev', 'splits L1 distances', 'splits max diff' );
  $plot1->gnuplot_plot_xy( $gens, #$topo_L1s, $topo_max_diff, 
			   $splits_avg_stddevs, $splits_L1s, $splits_max_diff);


  # inverse effective sample size
  $plot1->gnuplot_cmd(" set ytics ( $no_ytics_labels_string ) ");
  $plot1->gnuplot_cmd(" set origin $scale, $scale ");
  $plot1->gnuplot_set_plot_titles( 'tree length 1/ESS',
				   'alpha 1/ESS', 'pinvar 1/ESS' );
  $plot1->gnuplot_plot_xy( $gens, \@TL_invESSs, \@alpha_invESSs,
			   \@pinvar_invESSs );


  # lower two plots

  # Kolmogorov-Smirnov D statistic
  $plot1->gnuplot_cmd(" set xtics ( $xtics_labels_string ) ");
  $plot1->gnuplot_cmd(" set yrange [0.01:1] ");
  $plot1->gnuplot_cmd(" set ytics auto "); # $ytics_string ");

  $plot1->gnuplot_cmd(" set origin 0,0 ");
  $plot1->gnuplot_set_plot_titles( 'LogL KSD', 'tree length KSD',
				   'alpha KSD', 'p_{invar} KSD' );
  $plot1->gnuplot_plot_xy( $gens, $LnL_KSDs, $TL_KSDs, $alpha_KSDs,
			   $pinvar_KSDs );

  # Avg L1 distance
  $plot1->gnuplot_cmd(" set ytics ( $no_ytics_labels_string ) ")
    ;		    # ('x' 0.1 0, 'y' 0.01 0) " ); # $ytics_string ");
  $plot1->gnuplot_cmd(" set origin $scale, 0.0 ");
  $plot1->gnuplot_set_plot_titles(
				  'LogL avg L1',
				  'tree length avg L1',
				  'alpha avg L1',
				  'p_{invar} avg L1'
				 );
  $plot1->gnuplot_plot_xy( $gens, $LnL_L1s, $TL_L1s, $alpha_L1s,
			   $pinvar_L1s );

  $plot1->gnuplot_cmd(' unset multiplot ');
  # histograms



# topology histograms
  my $plot2 = Graphics::GnuplotIF->new( persist => $persist );
  $plot2->gnuplot_cmd(' set terminal x11 enhanced size 640,480 ') if ($enhanced);
  #    $plot2->gnuplot_cmd(' set multiplot ');
  $plot2->gnuplot_cmd(' set style data histeps ');
  my $histogram_filename = $base_filename . ".topology_histograms";
  my $histo_command = ' plot [-0.5:50.5] "' . $histogram_filename . '" using 2 ' . 't"run 1", ';
  for my $i_run (2..$n_runs) {
    $histo_command .= ' "" ' . ' using ' . ($i_run+1) . ' t"run ' . $i_run .'", ';
  }
  $histo_command .= ' "" using ($' . ($n_runs+2) . "/$n_runs" . ') t"average" ';
  $plot2->gnuplot_cmd( $histo_command );


# splits histograms
 my $plot2a = Graphics::GnuplotIF->new( persist => $persist );
  $plot2a->gnuplot_cmd(' set terminal x11 enhanced size 640,480 ') if ($enhanced);
  #    $plot2a->gnuplot_cmd(' set multiplot ');
  $plot2a->gnuplot_cmd(' set style data histeps ');
  $histogram_filename = $base_filename . ".splits_histograms";
  $histo_command = ' plot [-0.5:50.5] "' . $histogram_filename . '" using 2 ' . 't"run 1", ';
  for my $i_run (2..$n_runs) {
    $histo_command .= ' "" ' . ' using ' . ($i_run+1) . ' t"run ' . $i_run .'", ';
  }
  $histo_command .= ' "" using ($' . ($n_runs+2) . "/$n_runs" . ') t"average" ';
  $plot2a->gnuplot_cmd( $histo_command );

# histograms for LnL, TL, alpha, pinvar
  my $plot3 = Graphics::GnuplotIF->new( persist => $persist );
  $plot3->gnuplot_cmd(' set terminal x11 enhanced size 800,600')
    if ($enhanced);
  $plot3->gnuplot_cmd(' set style data histeps ');
  $plot3->gnuplot_cmd(' set multiplot ');
  $plot3->gnuplot_cmd(" set size 0.5,0.5 ");

  $histo_command = ' using 1:2 ' . 't"run 1", ';
  for my $i_run (2..$n_runs) {
    $histo_command .= ' "" ' . ' using 1:' . ($i_run+1) . ' t"run ' . $i_run .'", ';
  }
  $histo_command .= ' "" using 1:($' . ($n_runs+2) . "/$n_runs" . ') t"average" ';

  # LnL
  $plot3->gnuplot_cmd(" set origin 0.0,0.5 ");
  $histogram_filename = $base_filename . ".LnL_histograms";
  $plot3->gnuplot_set_xlabel(' log likelihood ');
  $plot3->gnuplot_cmd( ' plot "' . $histogram_filename . '" ' . $histo_command);

  # TL
  $plot3->gnuplot_cmd(" set origin 0.5,0.5 ");
  $histogram_filename = $base_filename . ".TL_histograms";
  $plot3->gnuplot_set_xlabel(' tree length ');
 $plot3->gnuplot_cmd( ' plot "' . $histogram_filename . '" ' . $histo_command);
 
  # alpha
  $plot3->gnuplot_cmd(" set origin 0.0,0.0 ");
  $histogram_filename = $base_filename . ".alpha_histograms";
  $plot3->gnuplot_set_xlabel(' alpha ');
  $plot3->gnuplot_set_plot_titles( 'run 1', 'run 2', 'run 3', 'total' );
 $plot3->gnuplot_cmd( ' plot "' . $histogram_filename . '" ' . $histo_command);

  # pinvar
  $plot3->gnuplot_cmd(" set origin 0.5,0.0 ");
  $histogram_filename = $base_filename . ".pinvar_histograms";
  $plot3->gnuplot_set_xlabel(' p_{invar} ');
  $plot3->gnuplot_set_plot_titles( 'run 1', 'run 2', 'run 3', 'total' );
 $plot3->gnuplot_cmd( ' plot "' . $histogram_filename . '" ' . $histo_command);
 
  $plot3->gnuplot_cmd(" unset multiplot ");

  if($n_temps <= 1){
     $plot3->gnuplot_pause(0);
     exit;
}
  ######################## temperature swapping acceptance rates: ##########################
  my $mc3swap_filename = $base_filename . ".mc3swap";

  my $n_cols_per_run = int( $n_temps * ( $n_temps - 1 ) / 2 );
  my $col            = 2;
  my $w_wind         = $n_runs * 300;
  my $ylabel_space = 0.032;
  my $w_plot         = (1-$ylabel_space) / $n_runs;

  my $plot4 = Graphics::GnuplotIF->new( persist => $persist, style => 'points' );
  $plot4->gnuplot_cmd(" set terminal x11 enhanced size $w_wind,360 ")
    if ($enhanced);
  $plot4->gnuplot_cmd(" set log ");
  $plot4->gnuplot_cmd(" set multiplot ");
 
  my $origin_x = 0;
$plot4->gnuplot_set_ylabel("mc3 swap acceptance rate");
  for my $i_run ( 1 .. $n_runs ) {
    $plot4->gnuplot_cmd(" set origin $origin_x,0.0 ");
my $plot_width = $w_plot + (($i_run==1)? $ylabel_space : 0);
 $plot4->gnuplot_cmd(" set size $plot_width,1 ");
    my $plot_cmd       = "plot [][0.01:1] ";
    my $run_col_offset = 1 + $n_cols_per_run * $i_run;
    for my $i_gap ( 1 .. $n_temps - 1 ) {
      for my $i_lo_temp ( 1 .. $n_temps - $i_gap ) {
	my $ptitle =
	  "T" . ( $i_lo_temp - 1 ) . ":T" . ( $i_lo_temp - 1 + $i_gap );
	$plot_cmd .=
	  " '$mc3swap_filename' using " . '1:($' 
	    . $col . "/" . '$'
              . ( $col + 1 ) . ') t"'
		. $ptitle . '", ';
	$col += 2;
	$plot4->gnuplot_cmd($plot_cmd);
      }
    }				# loop over T gaps
    $plot_cmd =~ s/,\s*$//;	# remove final comma
    $plot4->gnuplot_cmd($plot_cmd);
    $origin_x += ($i_run ==1)? $w_plot+$ylabel_space : $w_plot;
    $plot4->gnuplot_cmd(" unset ylabel");
  }				# loop over runs
  $plot4->gnuplot_cmd(" unset multiplot ");
  $plot4->gnuplot_pause(0);
}
