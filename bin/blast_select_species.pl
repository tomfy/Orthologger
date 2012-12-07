#!/usr/bin/perl -w
use strict;

# get all the lines of a blast m8, or mcl abc file
# with (at least) one sequence being in a particular species.
# (both file formats have seq ids in the first two columns)

# usage example:
# blast_select_species.pl 'Solanum' 21species.gg < 21species.abc
# This will actually get everything with Solanum (e.g. lycopersicum, tuberosum, ...)

my $species = shift;

#my $blastm8file = shift;
my $ggfile = shift;

die "2 required arguments missing, or can't open file.\n" 
. "Usage example: blast_select_species 'Solanum_lycopersicum' 21species.gg < 21species.abc > Solyc.abc \n" if(!defined $species  or ! -f $ggfile);   

my %id_species = ();

my $grepout = `grep $species $ggfile`;
my @greplines = split("\n", $grepout);
foreach my $line (@greplines){
my @ids = split( " ", $line );
my $ggspecies = shift @ids;
$ggspecies =~ s/:$/ /;
warn "$species  $ggspecies dont agree.\n";
for (@ids) {
    $id_species{$_} = $species;
}
}

while (<>) {
    my @cols = split( " ", $_ );
    if ( exists $id_species{ $cols[0] } or exists $id_species{ $cols[1] } ) {
        print;
    }
}
