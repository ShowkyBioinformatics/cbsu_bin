#!/usr/bin/perl

use strict;
use warnings;

my $input_file = $ARGV[0];

open (FILE1, $input_file) or die ("Could not open file \n");

while ($input_file = <FILE1>){
chomp $input_file;
my @file_bits = split(/\s+/,$input_file);
my $gene = $file_bits[0];
my $NS_POLY = $file_bits[1];
my $NS_FIX = $file_bits[2];
my $S_POLY = $file_bits[3];
my $S_FIX = $file_bits[4];
my $codons = $file_bits[5];
my $FET = $file_bits[6];



#$PNPS = ($NS_POLY/$S_POLY);



my $PNPS = $S_POLY ? $NS_POLY / $S_POLY : 0;

#$DNDS = ($NS_FIX/$S_FIX);
	
my $DNDS = $S_FIX ? $NS_FIX/$S_FIX :0;

#my $NI = ($PNPS/$DNDS);
my $NI = $DNDS ? $PNPS/$DNDS :0;

my $final_NI = sprintf("%.4f",$NI);
	

my $alpha = (1-$final_NI);

my $final_FET = sprintf("%.4f",$FET);

print "$gene\t$NS_POLY\t$S_POLY\t$NS_FIX\t$S_FIX\t$codons\t$final_NI\t$alpha\t$final_FET \n";

	}

exit;


