#!/usr/bin/perl -w
use strict;

my $file_basename = shift;

# my $filename = $file_basename . ".trprobs";

print put_ids_back_in_newick($file_basename);

exit;



sub put_ids_back_in_newick{
	my $file_basename = shift;
	my $filename = $file_basename . ".trprobs";
# print $filename, "\n";
	my $ready = 0; 
	my %number_id = ();
	my $count_leaves = 0;
	my $count_trees = 0;
	my $post_prob_sum = 0;

	open my $fh, "<$filename";
	my $out_string = '';
	while(<$fh>){
		if(/begin trees/){
			$ready++;
#			print "found begin trees.\n";
		}elsif(/^\s+translate\s*$/){
			$ready++;
#			print "found translate. ready: $ready \n";
		}
		next if($ready < 2);
		if(/\s*(\d+)\s+(\S+)[,;]/){
#			print "seq number, id: $1 $2\n";
			$number_id{$1} = $2; 
			$count_leaves++;
		}elsif(/\s*tree tree_(\d+).*[&]W\s+(\d+[.]\d+)\]\s+([0-9,()]+);/){
			my $tree_topo_number = $1;
			$count_trees++;
			my $tree_topo_post_prob = $2;
			$post_prob_sum += $tree_topo_post_prob;
			my $tree_newick = $3;
			foreach my $seq_number (keys %number_id){
				my $seq_id = $number_id{$seq_number};
				$tree_newick =~ s/([(,])$seq_number([,)])/$1$seq_id$2/;
			}

#		print"$count_leaves $count_trees $post_prob_sum\n";
		$out_string .=  "$tree_topo_number  $tree_topo_post_prob newick  $tree_newick\n";
		}
	}
	close $fh;
	return $out_string;
}
