package Grouper;
use Scalar::Util qw (blessed);
use Hash::Ordered;
use List::Util qw (min max sum);
use TomfyMisc qw (median);
use 5.012;    # so can use each with arrays

=pod
each sequence belongs to a 'species', each species belongs to a group
construct from file specifying groups, for example:

groupname asterids
    Solanum_lycopersicum
    Coffea_arabica
groupname rosids
    Arabidopsis_thaliana
    Populus_trichocarpa
groupname disallowed
    Oryza_sativa
    Amborella_trichopoda

first from this construct a Hash::Ordered of Hash::Ordered
with keys being groupnames (in same order as in file), and values
being Hash::Ordered's with keys being species names (in same order as in file)
this is then accessible with get_groupnames_species 
Also construct a regular (not ordered) hash species_groupname giving the
groupname for each species and accessible with get_species_groupname.

=cut

sub new{
   my $class = shift;
   my $arg1 = shift;            # a file with groups info in it
   my $args = {};
   my $self  = bless $args, $class;

   my $groupname_species; # a Hash::Ordered of Hash::Ordered
   ######  get the groups of species  ######
   if (-f $arg1) {             # get the groups of species from a file
      my $groups_filename = $arg1;
      open my $fh_gsp, "<", "$groups_filename";
      $groupname_species = Hash::Ordered->new();
      my $groupname = undef;
      while (my $line = <$fh_gsp>) {
         next if($line =~ /^\s*#/); # skip lines which have only a comment
         if ($line =~ /^\s*groupname\s+(\S+)/) { # start of a new group
            $groupname = $1;
            $groupname_species->set($groupname => Hash::Ordered->new());
         } elsif (defined $groupname) { # add a species to the group
            $line =~ s/#.*$//;
            $line =~ s/,/ /g;
            for (split(" ", $line)) {
               $groupname_species->get($groupname)->set($_ => ''); # keys: species, values: '' 
            }
         }
      }
   } elsif ( blessed $arg1  and  $arg1->isa('Hash::Ordered')) {
      #   print STDERR "Constructing Grouper obj from Hash::Ordered of Hash::Ordered \n";
      $groupname_species = $arg1;
   } else {
      die "arg $arg1 to Grouper::new is neither filename nor Hash::Ordered.\n";
   }
   $self->initialize($groupname_species);
   return $self;
}

sub initialize{
   my $self = shift;
   my $groupname_species = shift;
   $self->{groupname_species} = $groupname_species; # same as $category_species
   $self->{species_groupname} = {}; # ref to plain hash (not ordered)
   $self->{grpname_spcount} = Hash::Ordered->new();
   $self->{grp_summary_info} = {};
   if ($groupname_species->isa('Hash::Ordered')) {
      for my $grp_name ($self->{groupname_species}->keys()) {
         $self->{grpname_spcount}->set($grp_name => Hash::Ordered->new());
         for my $sp ($self->{groupname_species}->get($grp_name)->keys()) {
            my $spgroups = $self->{species_groupname}->{$sp} //  [];
            push @$spgroups, $grp_name;
            $self->{species_groupname}->{$sp} = $spgroups;
            $self->{grpname_spcount}->get($grp_name)->set($sp => 0);
         }
      }
   }
}

sub get_group_species_string{
   my $self = shift;
   my $gss = '';
   for my $grpname ($self->{groupname_species}->keys()) {
      my $species = $self->{groupname_species}->get($grpname);
      for my $sp ($species->keys()) {
         $gss .= sprintf("# %20s  %28s \n", $grpname, $sp);
      }
   }
   return $gss;
}

#
sub group_species_venn{
   my $self = shift;
   my $group_name = shift;
   my $species1 = shift;
   my $species2 = shift;
   my $sp_vennregion = {};      # keys: species, values A,B,AB

# print join(", ", @$species1), "\n";
# print join(", ", @$species2), "\n";
   my %AorBspecies_count = ();
   my %Aspecies_count = ();
   for my $sp (@$species1) {
      if ($self->get_groupname_species()->get($group_name)->exists($sp)) {
         $AorBspecies_count{$sp}++;
         $sp_vennregion->{$sp} = 'A';
         $Aspecies_count{$sp}++;
      }
   }
   my %AandBspecies_count = ();
 my %Bspecies_count = ();
   for my $sp (@$species2) {
      if ($self->get_groupname_species()->get($group_name)->exists($sp)) {
         $AorBspecies_count{$sp}++;
         $Bspecies_count{$sp}++;
         if (exists $sp_vennregion->{$sp} and $sp_vennregion->{$sp} eq 'A') {
            $sp_vennregion->{$sp} = 'AB';
            $AandBspecies_count{$sp}++;
         } else {
            $sp_vennregion->{$sp} = 'B';
         }
      }
   }
   while (my($k,$v) = each %$sp_vennregion) {
      #print "$k  $v \n";
   }
   my $AandB_count = scalar keys %AandBspecies_count;
   my $AxorB_count = (scalar keys %AorBspecies_count) - $AandB_count;
   return (scalar keys %Aspecies_count, scalar keys %Bspecies_count, $AandB_count, $AxorB_count, $sp_vennregion);
}



sub group_species_counts{
   my $self = shift;
   my $arg1 = shift;
 #  my $id_species = shift;
   my $species_count = {};      # plain (not ordered) hash.
   if(ref $arg1 eq 'ARRAY'){ # array ref of species (e.g. from get_implicit_species
  #    print join(" ", @$arg1), "\n";
      for my $sp (@$arg1){
         $species_count->{$sp}++;
      }
      # while(my ($k, $v) = each %$species_count){
      #    print "$k  $v \n";
      # }
   }elsif (ref $arg1 eq 'HASH') { # $arg1 is hashref of  species:count pairs
         $species_count = $arg1;
   }else {
      die "arg1 should be hashref, or array ref of species names, but is ", ref($arg1), "\n";
   }
   # all values of %species_count will be integers > 0

   my $grp_sp_count = Hash::Ordered->new();
   for my $grp_name ($self->{groupname_species}->keys()) { # looks like should still work if sp belongs to multiple groups.
      my $this_grp_n_sp_present = 0; # counts number of species in this group which are present.
      my $this_grp_sp_count = Hash::Ordered->new(); # counts number of sequences of each species in this group.
      for my $sp ($self->{groupname_species}->get($grp_name)->keys()) {
         $this_grp_sp_count->set($sp => $species_count->{$sp} // 0);
         #      print "sp: $sp ", $species_count{$sp}, "\n";
         $this_grp_n_sp_present++ if(($species_count->{$sp} // 0) > 0);
      }
  #    print "BBB:  $grp_name    ",  join(",", $this_grp_sp_count->as_list()), "\n";
      $self->{grpname_spcount}->set($grp_name => $this_grp_sp_count);
      my @sp_seq_counts_in_group = $self->{grpname_spcount}->get($grp_name)->values();
      #print "CCC:   $grp_name   ", join(",", @sp_seq_counts_in_group), "\n";
      my ($med_nseq, $max_nseq, $total_nseq) = (
                                                median(@sp_seq_counts_in_group),
                                                max(@sp_seq_counts_in_group), 
                                                sum(@sp_seq_counts_in_group),
                                               );

  #   print "DDD:  $grp_name   ", join("; ", ($this_grp_n_sp_present, $med_nseq, $max_nseq, $total_nseq)), "\n";
      $self->{grp_summary_info}->{$grp_name} = [ $this_grp_n_sp_present, $med_nseq, $max_nseq, $total_nseq ];
   }
}

sub get_n_species_in_group{
   my $self = shift;
   my $groupname = shift;
   return scalar $self->get_group($groupname)->keys();
}


sub get_groupname_species{
   my $self = shift;
   return $self->{groupname_species};
}

sub group_speciescount_strings{
   my $self = shift;
   my $group_name = shift;

   my $sp_count = $self->{grpname_spcount}->get($group_name); # Hash::Ordered

   my $string = '';
   for ($sp_count->values()) {
      $string .= sprintf("%1i ", $_);
   }
   $string =~ s/\s+$//;
   my $summary_string = sprintf("%1i %2.1f %1i %1i", @{$self->{grp_summary_info}->{$group_name}});
   return ($summary_string, $string);
}

sub get_group_summary_info{
   my $self = shift;
   my $group_name = shift;
   #   print "XXXX: $group_name \n";
   return @{$self->{grp_summary_info}->{$group_name}};
}

sub get_group_sequences_string{
   my $self = shift;
   my $node = shift;
   my $group_name = shift;
   my $id_species = shift;

   my $ids = $node->get_implicit_names();
   my @grpids = ();
   for my $id (@$ids) {
      my $sp = $id_species->{$id};
      if ($self->get_groupname_species()->get($group_name)->exists($sp)) {
         push @grpids, $id;
      }
   }
   return join(",", @grpids);
}

# sub get_groupname_species{
#    my $self = shift;
#    return $self->{groupname_species};
# }

sub get_species_groupname{
   my $self = shift;
   return $self->{species_groupname};
}

# sub get_species_group{ # get the 
#    my $self = shift;
#    return $self->get_group($self->{species_groupname});
# }

sub get_groupnames{
   my $self = shift;
   return $self->{groupname_species}->keys();
}

sub get_ith_groupname{
   my $self = shift;
   my $index = shift // 0;
   my @names = $self->{groupname_species}->keys();
   #   print STDERR "group names: ", join(", ", @names), "\n";
   my $grp_name = $names[$index];
   return $grp_name;            # ($grp_name, $grp);
}

sub get_group{ # for the group with name $groupname, return the value, a Hash::Ordered containing species names (as keys)
   my $self = shift;
   my $groupname = shift;
   my $group = undef;
  #   print STDERR "Group name: $groupname \n";
   if (defined $groupname) {
      if (defined $self->{groupname_species}->get($groupname)) {
         $group = $self->{groupname_species}->get($groupname);
      } else {
         die "No group for groupname $groupname.\n";
      }
   } else {              # if no groupname argument, return 1st group.
      my @group_names = $self->{group}->values();
      $group = $group_names[0];
   }
   return $group;
}

sub species_and_group_counts{
   # for the sequence ids specified in the array ref argument $sequenceids:
   # number of leaves of each species ($species_leafcount)
   # number of species present of each group ($groupname_speciescount)
   # number of leaves present from each group ($groupname_leafcount)
   # hash of ids present ($subtree_id_presence). keys are ids in subtree, values are all 1.
   my $self = shift;
   my $sequenceids = shift; # array ref
   my $id_species = shift;
   my $species_groupname = $self->get_species_groupname(); # hashref; keys: species names (e.g. 'Coffea_arabica'), values: groupnames (e.g. 'asterids')
   my $species_leafcount = {}; # keys: species names, values: number of leaves of that species in subtree.
   my $groupname_speciescount = {}; # keys: groupnames, values: numbers of species in group which are present in subtree.
   my $groupname_leafcount = {};
   my $subtree_id_presence = {}; # keys: ids in subtree, values: 1

   for my $an_id (@$sequenceids){ # @{$node->get_implicit_names()}) {
      $subtree_id_presence->{$an_id} = 1;
      my $species = (exists $id_species->{$an_id})? $id_species->{$an_id} : die "id: [$an_id]; species is unknown.\n";
      my $groupnames = $species_groupname->{$species} // ['other'];
      for my $groupname (@$groupnames){
         $groupname_leafcount->{$groupname}++;
         $groupname_speciescount->{$groupname}++ if(! exists $species_leafcount->{$species});
      }
      $species_leafcount->{$species}++;
   }
   return ($species_leafcount, $groupname_speciescount, $groupname_leafcount, $subtree_id_presence);
}

sub ids_present_string{
   my $self = shift;
   my $seqid_species = shift;
   my $seqid_present = shift;   # href, keys: sequence ids, values: 1
   my $group_name = shift;      # only do this group
   print "group name: $group_name \n";
   my @species_str = $self->{groupname_species}->get($group_name)->as_list();
   my $species_idsstring = Hash::Ordered->new(@species_str);
   # my %groupname_idsstring = ();
   my @sorted_ids = sort keys %$seqid_present;
   for my $an_id (@sorted_ids) {
      my $sp = $seqid_species->{$an_id};
      if ($self->get_groupname_species()->get($group_name)->exists($sp)) {
         $species_idsstring->concat($sp, "$an_id,");
      }
   }
   my $outstring = "$group_name" . '[';
   for my $sp ($species_idsstring->keys()) {
      my $s = $species_idsstring->get($sp); # // '';
      $outstring .= "$sp:" . $s if($s);
      $outstring =~ s/,\s*$/;/;
   }
   $outstring .= ']';
   return $outstring;
}

sub three_subtrees_species_counts{ # doesn't belong in Grouper
   my $self = shift;
   my $node = shift;
   my $sequenceid_species = shift;
   my @children = $node->get_children();
   my $nchild = scalar @children;
   if ($nchild == 2) {
      my ($Lchild, $Rchild) = @children[0,1]; #  L subtree root, is $children[0], R is $children[1]
      my ($species_countL, $cat_spcountL, $cat_leafcountL, $id_presenceL) = $self->species_and_group_counts($Lchild->get_implicit_names(), $sequenceid_species);
      my ($species_countR, $cat_spcountR, $cat_leafcountR, $id_presenceR) = $self->species_and_group_counts($Rchild->get_implicit_names(), $sequenceid_species);
      my ($species_countA, $cat_spcountA, $cat_leafcountA, $id_presenceA) = $self->species_and_group_counts($node->get_tree()->get_root()->get_implicit_names(), $sequenceid_species);

      while (my($sp,$c) = each %$species_countL) {
         $species_countA->{$sp} -= $c;
      }
      while (my($sp,$c) = each %$species_countR) {
         $species_countA->{$sp} -= $c;
      }

while (my($sp,$c) = each %$species_countA) {
         if ($c < 0) {
            warn "$sp has count of $c in parent subtree.\n";
         }
      }
      $cat_spcountA = {};
      for my $sp (keys %$species_countA){
         my $count = $species_countA->{$sp} // 0;
         my $groupnames = $self->get_species_groupname()->{$sp} // ['other'];
         for my $groupname (@$groupnames){
    #     print "Species, group: $sp  $groupname  $count \n";
          $cat_spcountA->{$groupname}++ if($count > 0);
       }
      }
      return ($species_countL, $cat_spcountL, $cat_leafcountL, $id_presenceL, 
              $species_countR, $cat_spcountR, $cat_leafcountR, $id_presenceR, 
              $species_countA, $cat_spcountA, $cat_leafcountA, $id_presenceA);

   } elsif ($nchild == 0) {
      # node is leaf
   } elsif ($nchild > 2) {
      die "Tree is not binary.\n";
   } elsif ($nchild == 1) {
      warn "Node has exactly 1 child.\n";
   } else {
      die "Node has unexpected number of children: $nchild \n";
   }

}

1;
