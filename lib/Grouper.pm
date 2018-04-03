package Grouper;
use Scalar::Util qw (blessed);
use Hash::Ordered;
use List::Util qw (min max sum);
use TomfyMisc qw (median);
use 5.012;                      # so can use each with arrays

# each sequence belongs to a 'species', each species belongs to a group
# we start with a set of groups, specified by a Hash::Ordered of Hash::Ordered
# e.g. could be Hash::Ordered->new('group1' => Hash::Ordered('Coffea_arabica' => 1, 'Coffea_canephora' => 1), 'group2' => Hash::Ordered('Vitis_vin' => 1, 'Oryza_sat' => 1))
# 'new' constructs from this a hash whose keys are species, values are groups ( i.e. group names such as 'group1'  )


sub new {
   my $class = shift;
   my $groups = shift; # Hash::Ordered; keys are group names ('group1', ...); values are Hash::Ordereds of species in group (values are genome sizes (?))
   my $args = {};
   my $self  = bless $args, $class;
   $self->{groups} = $groups;
   $self->{species_groups} = {};
   $self->{grp_spcount} = Hash::Ordered->new();
   $self->{grp_summary_info} = {};
   # print "ref groups: ", ref $groups, "   ", $groups->isa('Hash::Ordered'), "\n";
   if ($groups->isa('Hash::Ordered')) {
  #    print "ref groups: ", ref $groups, "   ", $groups->isa('Hash::Ordered'), "\n";
      for my $grp_name ($self->{groups}->keys()) {
   #      print "$grp_name \n";
         $self->{grp_spcount}->set($grp_name => Hash::Ordered->new());
         my $grp = $self->{groups}->get($grp_name);
         for my $sp ($grp->keys()) {
    #        print "$sp ";
            $self->{species_groups}->{$sp} = $grp_name;
            $self->{grp_spcount}->get($grp_name)->set($sp => 0);
         }
  #       print "\n";
      }
   }
   return $self;
}


sub group_species_leafcounts{
   my $self = shift;
   my $node = shift;

   my %species_count = ();

   for my $leaf_species (@{$node->get_implicit_species()}) {
      $species_count{$leaf_species}++;
   }

   my $grp_sp_count = Hash::Ordered->new();
   for my $grp_name ($self->{groups}->keys()) {
      my $grp = $self->{groups}->get($grp_name);
      my $grp_sp_present = 0;
      for my $sp ($grp->keys()) {
         if (defined $species_count{$sp}) {
            $self->{grp_spcount}->get($grp_name)->set($sp => $species_count{$sp});
            $grp_sp_present++;
         }
      }
      my @sp_seq_counts_in_group = $self->{grp_spcount}->get($grp_name)->values();
#      print "  $grp_name  ", join(", ", @sp_seq_counts_in_group), "\n";
      my ($med_nseq, $max_nseq, $total_nseq) = (
                                                median(@sp_seq_counts_in_group),
                                                max(@sp_seq_counts_in_group), 
                                                sum(@sp_seq_counts_in_group),
                                               );
      $self->{grp_summary_info}->{$grp_name} = [ $grp_sp_present, $med_nseq, $max_nseq, $total_nseq ];
   }
}

sub group_speciescount_strings{
   my $self = shift;
   my $group_name = shift;

   my $sp_count = $self->{grp_spcount}->get($group_name);

   my $string = '';
   for ($sp_count->values()) {
      $string .= sprintf("%2i ", $_);
   }
   $string =~ s/\s+$//;
   my $summary_string = sprintf("%2i %3.1f %2i %2i", @{$self->{grp_summary_info}->{$group_name}});
   return ($summary_string, $string);
}

sub get_group_summary_info{
   my $self = shift;
   my $group_name = shift;
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
      my $grp_name = $self->{species_groups}->{$sp};
      push @grpids, $id if($grp_name eq $group_name);
   }
   return join(",", @grpids);
}


1;
