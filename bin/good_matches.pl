#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );

# 
#my $abc_file = shift;
my $seq_id_regex = shift || die "good_matches.pl called with no arguments. must supply sequence id regular expression as argument.\n"
 . "Usage example: good_matches.pl 'Medtr5g017860[.]1' 250 < Medtr.abc > Medtr5g017860.1.matches   \n"
. "Gets the best 250 matches to Medtr5g017860.1 in Medtr.abc .";
my $max_n_to_keep = shift || 1200;


#open my $fhin, "<$abc_file" || die "Could open $abc_file for reading. \n";

my %seqpair_eval = ();

while (<>) {
    my ( $id1, $id2, $eval ) = split( " ", $_ );

    my $id12_key;
    if ( $id1 =~ /$seq_id_regex/ ) {
        $id12_key = $id1 . "\t" . $id2;
    }
    elsif ( $id2 =~ /$seq_id_regex/ ) {
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
	printf("%-30s  %-30s  e: %g \n", split(" ", $_), $seqpair_eval{$_});
}
}
	
