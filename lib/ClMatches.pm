package ClMatches;
use strict;
use List::Util qw (min max sum);
use constant LOG10 => log(10.0);


# ClMatch object stores the qids and their matches
# then, when all have been stored, choose just the best matches

sub new {
   my $class = shift;
   my $qidstr = shift; # comma separated qids string
   my $min_sim = shift;
   my $max_sim = shift;
   my $args  = {clusterqids_str => $qidstr, 
                clusterqids_matchcount => {}, # key: qid, value: count of matches to qid
                id2_evs => {}, # key: id2, value: arrayref of evs
                min_sim => $min_sim, max_sim => $max_sim, max_ev => -1,
                };
   my $self  = bless $args, $class;
   my @qids = split(",", $qidstr);
   $self->{n_qids} = scalar @qids;
   for(@qids){
      $self->add_qid($_);
   }
   return $self;
}

sub add_qid{
   my $self = shift;
   my $qid = shift;
   $self->{clusterqids_matchcount}->{$qid} = 0;
}

sub add_match{
   my $self = shift;
   my $qid = shift;
   my $id2 = shift;
   my $ev = shift;
die "qid $qid \n" if(!defined $id2);
   $self->{id2_evs}->{$id2} = [] if(!exists $self->{id2_evs}->{$id2}); 
   $self->{max_ev} = max($ev, $self->{max_ev});
   push @{$self->{id2_evs}->{$id2}}, $ev;
   $self->{clusterqids_matchcount}->{$qid}++;
}

sub add_matches{
   my $self = shift;
   my $id2_ev_string = shift;
   my @id1id2evs = split("\n", $id2_ev_string);
#   print STDERR "ZZZZ: ", $id2_ev_string, "\n";
#   print STDERR "ZZZ, ", join(";", @id1id2evs), "\n";
   for my $id1id2ev (@id1id2evs){
      my ($qid, $id2, $ev) = split(" ", $id1id2ev);
      $self->add_match($qid, $id2, $ev)
   }
}

sub all_qids_done{ # returns 1 if matches have been stored for all qids
   my $self = shift;
   my $qids_done_count = 0;     # 
# print STDERR "\n", "top of qll_qids_done \n";
   while (my ($qid, $count) = each %{$self->{clusterqids_matchcount}}) {
  #    print "qid: $qid  count: $count \n";
      $qids_done_count++ if($count >= 1);
   }
#   print STDERR "$qids_done_count   ", $self->{n_qids}, "\n";
   my $done =  $qids_done_count == $self->{n_qids};
 #   print STDERR $self->{clusterqids_str}, "   ", $qids_done_count, "  ", $self->{n_qids}, "\n" if(!$done);
   return $done;
}

sub set_min_sim{
my $self = shift;
my $new_min_sim = shift;
$self->{min_sim} = $new_min_sim;
}

sub top_n_by_avg_sim{
   my $self = shift;
   my $n = shift;
   $self->set_min_sim(0.5*($self->ev_to_sim($self->{max_ev}) + $self->{min_sim}));
   my %id2_simavg = ();
   while (my ($id2, $evs) = each %{$self->{id2_evs}}) {
      my $nevals = scalar @$evs;
      for my $eval (@$evs) {
         my $sim = $self->ev_to_sim($eval);
         $id2_simavg{$id2} += $sim; # $self->ev_to_sim($eval);
      }
      $id2_simavg{$id2} = ($id2_simavg{$id2} + $self->{min_sim}*($self->{n_qids} - $nevals))/$self->{n_qids};
   }
   my @skeys = sort { ($id2_simavg{$b} <=> $id2_simavg{$a})  ||  ($a cmp $b) } keys %id2_simavg;
   my $top_n_matches_str = '';
   my $last = min($n, scalar @skeys) - 1;
   for (0..$last) {
      my $id2 = $skeys[$_];
      $top_n_matches_str .= sprintf("  %s %6.2f \n", $id2,  $id2_simavg{$id2});
   }
   return $top_n_matches_str;
}

sub ev_to_sim{
   my $self = shift;
   my $evalue = shift;
   my $sim = ($evalue == 0)? $self->{max_sim} : -log($evalue)/LOG10;
   return $sim;
}

sub get_qids_str{
   my $self = shift;
   return $self->{clusterqids_str};
}

sub clusterqids_matchcount_str{
   my $self = shift;
   my $str = '';
   while (my($k, $v) = each $self->{clusterqids_matchcount}) {
      $str .= "$k $v ;  ";
   }
   $str .= "\n";
   return $str;
}


1;
