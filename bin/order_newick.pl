#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );
my $n_taxa = shift;
my $selected_column = shift || 0;
my $mask = (1 << $n_taxa) - 1;
while (<>) {
my @cols = split(" ", $_);
  my $newick = $cols[$selected_column];;
chomp $newick;
#print "newick: $newick \n";
$newick =~ s/^[^(]*//;
$newick =~ s/[^)]*$//;
#print "newick: $newick \n";
  my %split_count = ();
  my ($minlabel, $onewick, $x, $y) = order_newick($newick, \%split_count, $mask, 0);
print "newick: $newick\n", "ordered newick: $onewick\n";

for (keys %split_count) {
  print "$_ ", $split_count{$_}, "\n";
}

}


sub order_newick {
  my $newick = shift;
  my $split_count = shift;
  my $mask = shift;
  my $depth = shift;
  if (!defined $depth) {
    $depth = 0;
  }
# print "$depth \n";
  #print STDERR "\n" if($depth == 0);
  #print STDERR "$depth $newick \n";
  exit if($depth > 10);
#  print "   newick: $newick \n";
  if
  # ( $newick =~ /^([^,]+)$/ ) { # no commas -> this is a leaf!
   ( $newick =~ /^(\d+)(:\d*[.]?\d+)?$/ ) { # this subtree is leaf!
#    print "leaf branch. newick: $newick \n";
    my $split_bp = 1 << ($1 - 1);
    return ( $1, $newick, $1, $split_bp );
  } else {			# subtree has > 1 leaf.
#    print "non-leaf branch. newick: $newick \n";
    my $split_bp = 0;
    my %label_newick = ();
    my %cladeleaflabelset_count = ();
    $newick =~ /^[(](.*)[)](:\d+[.]\d+)?$/;   # remove outer (), and final branch length: e.g. ':0.15'

    my @newick_chars = split( '', $1 );	# without surrounding ()
    my $lmr_paren_count = 0;	# count of left parens - right parens
    my ( $il, $ir ) = ( 0, 0 );
    my $n_chars   = scalar @newick_chars;
    my $min_label = 10000000;
    foreach (@newick_chars) {
      die "$_ ", $newick_chars[$ir], " not same!\n"
	if ( $_ ne $newick_chars[$ir] );
      if ( $_ eq '(' ) {	# encountered left paren.
	$lmr_paren_count++;
      }
      if ( $_ eq ')' ) {	# encountered right paren.
	$lmr_paren_count--;
      }


      if (   ( $ir == $n_chars - 1 )
	     or ( $_ eq ',' and $lmr_paren_count == 0 ) ) { #split
	my $ilast = ( $ir == $n_chars - 1 ) ? $ir : $ir - 1;
	my $sub_newick = join( '', @newick_chars[ $il .. $ilast ] );

#	       print "subnewick [$sub_newick] il, ilast: $il $ir \n";
	my ( $label, $ordered_subnewick, $clade_leaf_labels, $subtree_split_bp ) = order_newick( $sub_newick, $split_count, $mask, $depth + 1 );
	#print "AAA:  $label; [$ordered_subnewick]; $clade_leaf_labels; $subtree_split_bp \n";
	$label_newick{$label} = $ordered_subnewick;
	$cladeleaflabelset_count{$clade_leaf_labels}++;

	$split_bp |= $subtree_split_bp; # bitwise OR of all subtree bitpatterns
	$min_label = min( $min_label, $label );
	$il        = $ir + 1;
	$ir        = $il;	# skip the ','
#print "il, ir: $il $ir \n";
      } else {
	$ir++;
      }

    }				# loop over chars in @newick_chars
    my $ordered_newick = '';
    foreach ( sort { $a <=> $b } keys %label_newick ) {
      $ordered_newick .= $label_newick{$_} . ",";
    }
    $ordered_newick =~ s/,$//;
    $ordered_newick = '(' . $ordered_newick . ')';
    my $clade_leaf_labels = join("_", keys %cladeleaflabelset_count);

    my $split_bp_complement = ~$split_bp & $mask;
    my $std_split_bp = ($split_bp % 2 == 1)? $split_bp : ~$split_bp & $mask; ##_complement;
    $split_count->{$std_split_bp}++ if($std_split_bp != $mask);
    return ( $min_label, $ordered_newick, $clade_leaf_labels, $split_bp );
  }
  die "shouldnt get here, in order_newick\n";
}

