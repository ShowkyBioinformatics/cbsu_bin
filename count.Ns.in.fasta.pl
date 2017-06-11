#!/usr/local/bin/perl
use strict;
use Bio::SeqIO;

my $input_file = $ARGV[0];

open (FILE1, $input_file) or die ("Could not open file \n");

# edit fasta file name
my $seqio = Bio::SeqIO->new(-file => "$input_file", -format => 'fasta');
while(my $seq = $seqio->next_seq) {
  my $s = $seq->seq;
  $s =~ s/n/N/g;
  my $id  = $seq->display_id; # name (ID)
  my $len = $seq->length; # length
  my %count;
  $count{$_}++ foreach split //, $s; # count bases
  $count{"N"} = 0 if !exists($count{"N"}); # create N = 0 if no N in sequence
  my $percn = ($count{"N"}/$len)*100; # ratio of N to not N
  my @output = ($id, $len, $count{"N"}, $percn);
  print join("\t", @output), "\n";
}
