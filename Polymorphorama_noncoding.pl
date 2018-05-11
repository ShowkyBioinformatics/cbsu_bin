#!/usr/bin/env perl 

# requires complilation of C program - tajd.c
# gcc -o tajd tajd.c -lm

# produces common summaries of NONcoding data including the frequency spectrum of mutations
# to run: perl Polymorpho_NC.pl infile frequency-cut-off reconstr(Yes=1/No=0/use_list_names=-1)
# where the infile contains a list with names of NONA alignments to analyse
# reconstr(Yes=1/No=0) refers to whether the filename to be used contains "reconstr" - i.e. bc the ANC was reconstructed using PAML
# needs .nona formatted files, with one outgroup (listed first) followed by ingroup_poly
# also needs executable "tajd" to calculate Tajima's D and Fu & Li's D

my $infile = $ARGV[0];
my $freq_cut_off = $ARGV[1];
my $reconstr = $ARGV[2];

open (IN1, $infile);
open (OUT, '>frequencies.xls') or die "can't open outfile\n";
open (OUT2, '>summarystats.xls') or die "can't open outfile\n";

print OUT2 "locus \t outgroup \t samp_size \t NonCodingSites \t S_NC  \t S_NC_freq>$freq_cut_off \t D_NC  \t pi_NC  \t pi_JC_NC  \t Dxy_NC  \t Dxy_JC_NC \t Dxy_Kimura_NC  \t TajD_NC  \t Fu&LiD_NC \tFayWu_H \n";
my @filehandles = <IN1>;


# array of arrays from 1 to numseq with count of polymorphic variants in each frequency class (from 1 to numseq-1)
# i.e. a singleton is in frequency class $poly_freq_Syn[1]

foreach $file (@filehandles)
  {
    chomp $file;

    if ($reconstr ==1)
		{
		$filereconstr=$file.".reconstr";
		open IN2, "< $filereconstr" or die "wrong format for infile or file name must end with '.reconstr'!\n";
		}
    elsif ($reconstr ==0)
        {
		$file_1og = $file.".1og";
        open IN2, "< $file_1og" or die "wrong format for infile or file must en with 1og!\n";
        }
	else 
		{
	    open IN2, "< $file" or die "wrong format for infile - no ending for infile !\n";
		}	

	my @data = ();
	my @sequence_temp =();  my @data_temp =();	# needed to read in nexus files
	my @sequence_names = ();
    my @poly_freq_NC = ();	  # array from 1 to numseq with count of polymorphic variants in each frequency class (from 1 to numseq-1)
							  # i.e. a singleton is in frequency class $poly_freq_NC[1]
	my @freqNC_Ts = (); # non-coding transitions
	my @freqNC_Tv = (); # non-coding transversions

	$numseqs=0;

	# changed this part to do fasta alignments instead of nona
	LOOP: while (my $line = <IN2>)
		{
		last LOOP if $line =~ m/\;/go; #end the loop if the line contains a semicolon anywhere
		chomp ($line);
		$line =~ s/\s+//; # remove all white spaces
		if ($line =~ />/) 
			{
			push(@sequence_names, $line);
			$sequence_temporary = join ("", @data_temp);
			push(@sequence_temp,$sequence_temporary); 
			@data_temp=();
			} 
		else {push(@data_temp, $line);}
		}
	close IN2;
	$sequence_temporary = join ("", @data_temp);
	push(@sequence_temp,$sequence_temporary); 
	for ($x=0; $x< scalar @sequence_names; $x++)	# put sequences into @data
		{
		$data[$x]=$sequence_temp[$x+1];
		#print $data[$x], "\n";
		}


# print out original data
	# print "data were:\n";
	foreach (@data)
		{
		split(//, $data);
		# print $_,"\n";
		$numseqs++;
		}
	
	chomp $file;	
	print "\nlocus ", $file, " numseqs: ", $numseqs, "\n";
	print "outgroup used: ", $sequence_names[0], "\n";

#assign information to variable codons
	my $seqlen = length($data[0]);



# initialize arrays for polymorphisms ($numseqs+1 for divergence at position $numseqs+1)
for ($ind=0; $ind<($numseqs+1); $ind++)
   {
   $poly_freq_NC[$ind]=0;
   $poly_freq_NC_temp[$ind]=0;
   $freqNC_GC_AT[$ind]=0;
   $freqNC_AT_AT[$ind]=0;
   $freqNC_AT_GC[$ind]=0;
   $freqNC_GC_GC[$ind]=0;
   $freqNC_Ts[$ind]=0;
   $freqNC_Tv[$ind]=0;
   }

$no_NC_sites= $gc_content=0; $no_NC_sites_invariant=$gc_content_invariant=0;

# first task, let take polymorphic codons, and feed them to subrouting "codon_processor" if there is no gap
#**********************************************************************************************************************************************************
# big loop through all codons of a locus
for ($pos=0; $pos<$seqlen; $pos++)
	{
   
   # check for a gap  - only proceed if there is no gap
    $switch_gap=0;
    for ($i=0; $i<$numseqs; $i++)
       {
	   if ((substr($data[$i],$pos,1) eq ':') or (substr($data[$i],$pos,1) eq '-'))
		  {
		  $switch_gap=1;
		  }
	    }
		
  if ($switch_gap==0)		
	{
	#first count the number of noncoding sites if there is no gap
	$no_NC_sites++;

    for ($i=0; $i<$numseqs; $i++)
       {
	   if ((substr($data[$i],$pos,1) eq 'G') or (substr($data[$i],$pos,1) eq 'C'))
		  {
		  $gc_content++;
		  }
	    }

    # assing polymorphism and divergence
	@pos_array = ();
	for ($x=0; $x<$numseqs; $x++)
		{
		$pos_array[$x]=(substr($data[$x], $pos, 1));
		}
	$position = join ("", @pos_array);
	$position_OG= substr($position, 0, 1);
	
      #	print "\n";
	  #	print "position: ", $position, "\n";


	$freq_G=$freq_C=$freq_T=$freq_A=$freq_OG=0; 
	
	while ($position =~ /G/g){$freq_G++;}
	while ($position =~ /A/g){$freq_A++;}
	while ($position =~ /T/g){$freq_T++;}
	while ($position =~ /C/g){$freq_C++;}
	while ($position =~ /$position_OG/g){$freq_OG++;}

$most_common='G';
$freq_most_common=$freq_G;
if ($freq_A>$freq_most_common)
	{$most_common='A';	$freq_most_common=$freq_A;}
if ($freq_T>$freq_most_common)
	{$most_common='T';	$freq_most_common=$freq_T;}
if ($freq_C>$freq_most_common)
	{$most_common='C';	$freq_most_common=$freq_C;}

# print "G: $freq_G; C: $freq_C; A: $freq_A; T: $freq_T\n";

	if (($freq_G==$numseqs) or ($freq_C==$numseqs) or ($freq_T==$numseqs) or ($freq_A==$numseqs))  # monomorphic
		{
		# count GC-content for invariant sites only
		$no_NC_sites_invariant++;
		if (($freq_G==$numseqs) or ($freq_C==$numseqs)) {$gc_content_invariant++;}
		} 
	elsif (($freq_OG==1) and ($freq_most_common==($numseqs-1)))	#divergence
		{
		$poly_freq_NC[$numseqs]++;
		if (($position_OG eq 'G') or ($position_OG eq 'C'))	{if (($most_common eq 'G') or ($most_common eq 'C')) {$freqNC_GC_GC[$numseqs]++;} else {$freqNC_GC_AT[$numseqs]++;} }
		if (($position_OG eq 'A') or ($position_OG eq 'T'))	{if (($most_common eq 'A') or ($most_common eq 'T')) {$freqNC_AT_AT[$numseqs]++;} else {$freqNC_AT_GC[$numseqs]++;} }
		
		if ($position_OG eq 'A') {if ($most_common eq 'G') {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		if ($position_OG eq 'G') {if ($most_common eq 'A') {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		if ($position_OG eq 'C') {if ($most_common eq 'T') {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		if ($position_OG eq 'T') {if ($most_common eq 'C') {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		}
	elsif ($freq_OG>1) # polymorphic (no divergence)
		{
		if ($position_OG ne 'G') {$poly_freq_NC[$freq_G]++; if ($position_OG eq 'C') {$freqNC_GC_GC[$freq_G]++;} else {$freqNC_AT_GC[$freq_G]++;}}
		if ($position_OG ne 'A') {$poly_freq_NC[$freq_A]++; if ($position_OG eq 'T') {$freqNC_AT_AT[$freq_A]++;} else {$freqNC_GC_AT[$freq_A]++;}}
		if ($position_OG ne 'T') {$poly_freq_NC[$freq_T]++; if ($position_OG eq 'A') {$freqNC_AT_AT[$freq_T]++;} else {$freqNC_GC_AT[$freq_T]++;}}
		if ($position_OG ne 'C') {$poly_freq_NC[$freq_C]++; if ($position_OG eq 'G') {$freqNC_GC_GC[$freq_C]++;} else {$freqNC_AT_GC[$freq_C]++;}}
		
		if ($position_OG ne 'G') {if ($position_OG eq 'A') {$freqNC_Ts[$freq_G]++;} else {$freqNC_Tv[$freq_G]++;}}
		if ($position_OG ne 'A') {if ($position_OG eq 'G') {$freqNC_Ts[$freq_A]++;} else {$freqNC_Tv[$freq_A]++;}}
		if ($position_OG ne 'T') {if ($position_OG eq 'C') {$freqNC_Ts[$freq_T]++;} else {$freqNC_Tv[$freq_T]++;}}
		if ($position_OG ne 'C') {if ($position_OG eq 'T') {$freqNC_Ts[$freq_C]++;} else {$freqNC_Tv[$freq_C]++;}}
		}
	elsif (($freq_OG==1) and ($freq_most_common<($numseqs-1)))	#poly and divergence
		{
		$poly_freq_NC[$numseqs]++;
		if (($position_OG eq 'G') or ($position_OG eq 'C'))	{if (($most_common eq 'G') or ($most_common eq 'C')) {$freqNC_GC_GC[$numseqs]++;} else {$freqNC_AT_GC[$numseqs]++;} }
		if (($position_OG eq 'A') or ($position_OG eq 'T'))	{if (($most_common eq 'A') or ($most_common eq 'T')) {$freqNC_AT_AT[$numseqs]++;} else {$freqNC_GC_AT[$numseqs]++;} }
 
		if (($position_OG ne 'G') and ($most_common ne 'G')) {$poly_freq_NC[$freq_G]++; {if ($most_common eq 'C') {$freqNC_GC_GC[$freq_G]++;} else {$freqNC_AT_GC[$freq_G]++;} }}
		if (($position_OG ne 'A') and ($most_common ne 'A')) {$poly_freq_NC[$freq_A]++; {if ($most_common eq 'T') {$freqNC_AT_AT[$freq_A]++;} else {$freqNC_GC_AT[$freq_A]++;} }}
		if (($position_OG ne 'T') and ($most_common ne 'T')) {$poly_freq_NC[$freq_T]++; {if ($most_common eq 'A') {$freqNC_AT_AT[$freq_T]++;} else {$freqNC_GC_AT[$freq_T]++;} }}
		if (($position_OG ne 'C') and ($most_common ne 'C')) {$poly_freq_NC[$freq_C]++; {if ($most_common eq 'G') {$freqNC_GC_GC[$freq_C]++;} else {$freqNC_AT_GC[$freq_C]++;} }}
		
		if ($position_OG eq 'G')	{if ($most_common eq 'A')  {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		if ($position_OG eq 'A')	{if ($most_common eq 'G')  {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		if ($position_OG eq 'T')	{if ($most_common eq 'C')  {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		if ($position_OG eq 'C')	{if ($most_common eq 'T')  {$freqNC_Ts[$numseqs]++;} else {$freqNC_Tv[$numseqs]++;} }
		
		if (($position_OG ne 'G') and ($most_common ne 'G'))  {if ($most_common eq 'A') {$freqNC_Ts[$freq_G]++;} else {$freqNC_Tv[$freq_G]++;} }
		if (($position_OG ne 'A') and ($most_common ne 'A'))  {if ($most_common eq 'G') {$freqNC_Ts[$freq_A]++;} else {$freqNC_Tv[$freq_A]++;} }
		if (($position_OG ne 'T') and ($most_common ne 'T'))  {if ($most_common eq 'C') {$freqNC_Ts[$freq_T]++;} else {$freqNC_Tv[$freq_T]++;} }
		if (($position_OG ne 'C') and ($most_common ne 'C'))  {if ($most_common eq 'T') {$freqNC_Ts[$freq_C]++;} else {$freqNC_Tv[$freq_C]++;} }
		
		}
	else	
		{die "mutation not assigned\n";}
	     
	 
	# print "polytable NonCoding: ", join ("-", @poly_freq_NC), "\n";
	
	}  # loop for no_GAP --- if ($switch_gap==0)	
}  # loop for all codons $pos


$no_polyNC=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$no_polyNC=$no_polyNC+$poly_freq_NC[$ind];
	}
$no_divNC=$poly_freq_NC[$numseqs];
$no_div_GC_GC=$freqNC_GC_GC[$numseqs];	$no_div_AT_AT=$freqNC_AT_AT[$numseqs];	$no_div_AT_GC=$freqNC_AT_GC[$numseqs];	$no_div_GC_AT=$freqNC_GC_AT[$numseqs];	
$no_div_Ts=$freqNC_Ts[$numseqs]; $no_div_Tv=$freqNC_Tv[$numseqs];

if ($freq_cut_off != 0)
	{
	$no_polyNC_freq=0;
	for ($ind=1; $ind<$numseqs-1; $ind++)
		{
		if ( ($ind/($numseqs-1)) > $freq_cut_off)
			{
			$no_polyNC_freq=$no_polyNC_freq+$poly_freq_NC[$ind];
			}
		}
	}


# calculate pi
$pi_NC_total=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$pi_NC[$ind]=(2*($ind/($numseqs-1))*(1-($ind/($numseqs-1))))*$poly_freq_NC[$ind];
	$pi_NC_total=$pi_NC_total+$pi_NC[$ind];
	}
$pi_NC_total=$pi_NC_total*(($numseqs-1)/($numseqs-2));	

$pi_NC_site=$pi_NC_total/$no_NC_sites;

$pi_JC_NC= -0.75*log(1-(4/3)*$pi_NC_site);

# calculate Dxy
$freq_NC=0; $freq_NC_Ts=$freq_NC_Tv=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$freq_NC_temp[$ind]=$ind/($numseqs-1)*$poly_freq_NC[$ind];
	$freq_NC=$freq_NC+$freq_NC_temp[$ind];
	
	$freq_NC_Ts_temp[$ind]=$ind/($numseqs-1)*$freqNC_Ts[$ind];
	$freq_NC_Ts=$freq_NC_Ts+$freq_NC_Ts_temp[$ind];
	$freq_NC_Tv_temp[$ind]=$ind/($numseqs-1)*$freqNC_Tv[$ind];
	$freq_NC_Tv=$freq_NC_Tv+$freq_NC_Tv_temp[$ind];
	}
	
$Dxy_NC = ($no_divNC+$freq_NC)/$no_NC_sites;
$Dxy_NC_Ts = ($no_div_Ts+$freq_NC_Ts)/$no_NC_sites;
$Dxy_NC_Tv = ($no_div_Tv+$freq_NC_Tv)/$no_NC_sites;

$Dxy_JC_NC= -0.75*log(1-(4/3)*$Dxy_NC);
$Dxy_Kimura_NC=  0.5*log(1/(1 - (2*$Dxy_NC_Ts) - $Dxy_NC_Tv)) + (0.25*log(1/(1 - (2*$Dxy_NC_Tv))));

# calculate Fay&Wu's H
$pi_H_NC_total=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$pi_H_NC[$ind]=2*$ind*$ind*$poly_freq_NC[$ind];
	$pi_H_NC_total=$pi_H_NC_total+$pi_H_NC[$ind];
	}
$pi_H_NC_total=$pi_H_NC_total*(1/(($numseqs-1)*($numseqs-2)));
$FayWu_H=$pi_NC_total-$pi_H_NC_total;


#**************** calculate TajD

$numseqs_ingroup=$numseqs-1;
open INPUT, "./tajd ". $numseqs_ingroup . " " . $no_polyNC . " " . $pi_NC_total . " " .$poly_freq_NC[1] . " |";
# print "./tajd  $numseqs_ingroup  $no_polyNC  $pi_NC_total  $poly_freq_NC[1]\n";
	my $line1 = <INPUT>; chomp $line1;
	my $line2 = <INPUT>; chomp $line2;
	($junk, $TajD_NC)= split(/=/, $line1);
	($junk, $FuLiD_NC)= split(/=/, $line2);
	


pop @poly_freq_NC; pop @poly_freq_NC;
pop @freqNC_GC_GC; pop @freqNC_GC_GC; pop @freqNC_AT_AT; pop @freqNC_AT_AT; pop @freqNC_AT_GC; pop @freqNC_AT_GC; pop @freqNC_GC_AT; pop @freqNC_GC_AT;
pop @freqNC_Ts; pop @freqNC_Ts; pop @freqNC_Tv; pop @freqNC_Tv;
$poly_freq_NC[0]=$no_divNC;
$freqNC_GC_GC[0]=$no_div_GC_GC;	$freqNC_AT_AT[0]=$no_div_AT_AT; $freqNC_AT_GC[0]=$no_div_AT_GC; $freqNC_GC_AT[0]=$no_div_GC_AT;
$freqNC_Ts[0]=$no_div_Ts; $freqNC_Tv[0]=$no_div_Tv;

print OUT2 $file, "\t", $sequence_names[0],  "\t", $numseqs-1, "\t", $no_NC_sites,"\t", $no_polyNC, "\t", $no_polyNC_freq, "\t", $no_divNC, "\t",  $pi_NC_site, "\t", $pi_JC_NC, "\t", $Dxy_NC, "\t", $Dxy_JC_NC, "\t", $Dxy_Kimura_NC, "\t", $TajD_NC, "\t",  $FuLiD_NC, "\t", $FayWu_H,  "\n";

$gc_content=$gc_content/($no_NC_sites*($numseqs+1));
$gc_content_invariant=$gc_content_invariant/$no_NC_sites_invariant;

# print  "\nlocus ", $file, " numseqs: ", $numseqs, "\n";
# print "Num of NC Sites: ", $no_NC_sites, "\n";
# print  "NonCoding Poly: $no_polyNC, NonCoding Divergence: $no_divNC \n"; 
# print "Pairwise NonCoding diversity: $pi_NC_site \t JC corr:  $pi_JC_NC  \n";
# print "Dxy_NC ",  $Dxy_NC, "\tDxy_JC_NC  ", $Dxy_JC_NC , "\n";
# print "NonCoding Taj D is : $TajD_NC, Fu&Li D is $FuLiD_NC\n";
# print  "NonCoding Poly above a frequency of $freq_cut_off: $no_polyNC_freq  \n"; 
print OUT $file, "_NC\t" , join ("\t", @poly_freq_NC), "\n";
print OUT3 $file, "\t %GC(all_sites)\t" , $gc_content, "\n";
print OUT3 $file, "\t %GC(invariant_sites)\t" , $gc_content_invariant, "\n";
print OUT3 $file, "\t AT->GC\t" , join ("\t", @freqNC_AT_GC), "\n"; 
print OUT3 $file, "\t GC->AT\t" , join ("\t", @freqNC_GC_AT), "\n"; 
print OUT3 $file, "\t AT->AT\t" , join ("\t", @freqNC_AT_AT), "\n"; 
print OUT3 $file, "\t GC->GC\t" , join ("\t", @freqNC_GC_GC), "\n"; 

print OUT4 $file, "\t Ts\t" , join ("\t", @freqNC_Ts), "\n"; 
print OUT4 $file, "\t Tv\t" , join ("\t", @freqNC_Tv), "\n"; 


} # loop foreach file




