package FastaSelect;
use strict;
use warnings;

use Moose;
use namespace::autoclean;

has id_seq_hashref => (
		isa => 'HashRef',
		is => 'rw',
		default => sub{ {} }
		);

has fasta_filename => (
		isa => 'Str',
		is => 'ro',
		required => 1
		);


sub BUILD {
	my $self = shift;

	my $fastafilename = $self->fasta_filename();	
	open my $fasta_fh, "<$fastafilename";

	my %id_seq_hash = ();
	while(<$fasta_fh>){
		if(/^>/){
			chomp;
			my $idline = $_;
			my $sequence = <$fasta_fh>;
			chomp $sequence;

			if( $idline =~ /^>(\S+?)[\s]\s*(.*)/){
				my $id = $1;
				my $annotation = $2;
				$id_seq_hash{$id} = $2 . "\n" . $sequence;
			}
		}
	}
	close $fasta_fh;
	$self->id_seq_hashref(\%id_seq_hash);
}

sub select_fasta{
# argument $ids_to_select is either a filename
# of a file with ids to be selected in the first 
# (or zeroeth if you prefer) column, or
# a string with similar format.
# the id may be preceded in the line by whitespace
# or the > character immediately before the id.

	my $self = shift;
	my $ids_to_select = shift;
# print "[$ids_to_select]\n";
	my @ids = ();
	if( !($ids_to_select =~ /\n/) and  -f $ids_to_select){
		open my $ids_fh, "<$ids_to_select";

		while(<$ids_fh>){
			chomp;
			/^\s*>?(\S+)/; # exclude init whitespace and > if present, then take everything up to first whitespace.
				push @ids, $1;
		}
	}elsif(ref $ids_to_select eq 'HASH'){
		@ids = keys %$ids_to_select;
	}else{ # assume it is a string
		while($ids_to_select){
			$ids_to_select =~ s/^\s*>?(\S+).*[\n]?//;
			push @ids, $1;
		}
	}
	my $selected_fasta = '';
	my $ishr = $self->id_seq_hashref();
	foreach (@ids){
#	print "ID: $_\n";
		if( exists $ishr->{$_}){
#	print "found\n";
			my $annot_sequence = $ishr->{$_};
			my ($annot, $seq) = split("\n", $annot_sequence);
			$selected_fasta .= ">$_  $annot \n" . $seq . "\n";
		}else{
#	print "not found\n";
}
	}
	return $selected_fasta;
}

__PACKAGE__->meta->make_immutable;

1;
