#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );

# 
#my $abc_file = shift;
my $seq_id = shift || die "good_matches.pl called with no arguments. must supply sequence id regular expression as argument.\n"
. "Usage example: good_matches.pl 'Medtr5g017860[.]1' 250 < Medtr.abc > Medtr5g017860.1.matches   \n"
. "Gets the best 250 matches to Medtr5g017860.1 in Medtr.abc .";
my $max_n_to_keep = shift || 1200;

$seq_id =~ s/[.]/[.]/g;  # put .s in [] to show we only want to match the . character itself.

print STDERR "seq id regex to match: $seq_id \n";
#open my $fhin, "<$abc_file" || die "Could open $abc_file for reading. \n";

my %seqpair_eval = ();

while (<>) {
	my @cols = split(" ", $_);

	my ( $id1, $id2, $eval );

	if (scalar @cols == 12){
		($id1, $id2, $eval) = ($cols[0], $cols[1], $cols[10]);
	} elsif(scalar @cols == 3){
		($id1, $id2, $eval) =  ($cols[0], $cols[1], $cols[2]); 
	}
	else{
		die "File read from stdin does not appear to be either m8 (12 col) or abc (3 col) format.\n"
	}

	my $id12_key;
	if ( $id1 =~ /$seq_id/ ) {
		$id12_key = $id1 . "\t" . $id2;
	}
	elsif ( $id2 =~ /$seq_id/ ) {
		$id12_key = $id2 . "\t" . $id1;
	}
	else {
		next;    # match doesn't involve sequence of interest.
	}

	if ( exists $seqpair_eval{$id12_key} ) {
		my $id12_eval = $seqpair_eval{$id12_key};
		$id12_eval = min( $id12_eval, $eval );
		$seqpair_eval{$id12_key} = $id12_eval;
	}
	else {
		$seqpair_eval{$id12_key} = $eval;
	}

}

my @matches = sort {$seqpair_eval{$a} <=> $seqpair_eval{$b}} keys %seqpair_eval;

my $n_to_show = min($max_n_to_keep, scalar @matches);
if($n_to_show == 0){
	print STDERR "No matches found.\n";
}else{
	foreach (@matches[0..$n_to_show-1]){
		printf("%-30s  %-30s   %g \n", split(" ", $_), $seqpair_eval{$_});
	}
}

