#!/usr/bin/env perl                                                                                                                                                                                    
use warnings;
use strict;
use Bio::AlignIO;
use Bio::PopGen::IO;
use Bio::PopGen::Statistics;
use Bio::PopGen::Utilities;

my $infile = shift;

my $in = Bio::AlignIO->new(-format => 'fasta', -file => $infile);
my $aln;
if($aln = $in->next_aln){
}

my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment=>$aln, -include_monomorphic =>1);

my $pi = Bio::PopGen::Statistics->pi($pop);
my $theta = Bio::PopGen::Statistics->theta($pop);
my $D = Bio::PopGen::Statistics->tajima_D($pop);
my $segsites = Bio::PopGen::Statistics->segregating_sites_count($pop);
my $fu_and_li_D_star = Bio::PopGen::Statistics->fu_and_li_D_star($pop);
my $fu_and_li_F_star = Bio::PopGen::Statistics->fu_and_li_F_star($pop);

print "$pi"."\t"."$theta"."\t"."$D"."\t"."$fu_and_li_F_star"."\t"."$fu_and_li_D_star"."\t"."$segsites"."\n";
