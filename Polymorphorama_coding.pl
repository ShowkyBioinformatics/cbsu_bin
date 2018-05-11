#!/usr/bin/env perl 

# requires complilation of C program - tajd.c
# gcc -o tajd tajd.c -lm

# produces common summaries of CODING data including the frequency spectrum of mutations (synonymous, non-synonymous)
# performs codon-bias analysis using codon table for D. pseudoobscura

# to run: perl Polymorphorama1.pl infile frequency-cut-off reconstr_outgroup?(Yes=1/No=0/use_list_names=-1)
# where the infile contains a list with names of FASTA alignments to analyse
# 
# needs .nexus formatted files, with one outgroup (listed first) followed by ingroup_poly
# also needs executable "tajd" to calculate Tajima's D and Fu & Li's D

# frequency-cut-off ignores polymorphisms of a frequency of less than the cut-off for McDonald-Kreitman tables
# frequency-cut-off also implemented for codon bias stuff


my $infile = $ARGV[0];
my $freq_cut_off = $ARGV[1];

open (IN1, $infile);
system(qq{remcarm.pl $infile});

open (OUT, '>frequencies.xls') or die "can't open outfile\n";
open (OUT2, '>summarystats.xls') or die "can't open outfile\n";
open (OUT3, '>codon_bias.xls') or die "can't open outfile\n";

print OUT2 "locus \t outgroup \t samp_size \t SynSites \t S_Syn  \t S_Syn_freq>$freq_cut_off \t D_Syn \t pi_Syn  \t pi_JC_Syn  \t Dxy_Syn  \t Dxy_JC_Syn  \t D_Kimura_2para \t TajD_Syn  \t Fu&LiD_Syn \t Fay&WuH_Syn \t\t RepSites \t S_Rep  \t S_Rep_freq>$freq_cut_off \t D_Rep  \t pi_Rep  \t pi_JC_Rep  \t Dxy_Rep  \t Dxy_JC_Rep  \t TajD_Rep  \t Fu&LiD_Rep  \t Fay&WuH_Rep \n";
print OUT3 "locus \t\t S_freq>$freq_cut_off \t divergence  \t singleton  \t doubleton  \t tripleton  \n";

my @filehandles = <IN1>;

# array of arrays from 1 to numseq with count of polymorphic variants in each frequency class (from 1 to numseq-1)
# i.e. a singleton is in frequency class $poly_freq_Syn[1]

foreach $file (@filehandles)
  {
    chomp $file;
	open IN2, "< $file" or die "wrong format for infile !\n";
	
        
    my @data = (); 
	my @sequence_temp =();  my @data_temp =();	# needed to read in nexus files
    my @sequence_names = ();
    my @poly_freq_Syn = ();	  # array from 1 to numseq with count of polymorphic variants in each frequency class (from 1 to numseq-1)
							  # i.e. a singleton is in frequency class $poly_freq_Syn[1]
	my @poly_freq_Rep = ();	  # array from 1 to numseq with count of polymorphic variants in each frequency class (from 1 to numseq-1)

	my @freq_P_U = ();	# preferred to unpreferred
	my @freq_U_P = ();	# unpreferred to preferred
	my @freq_P_P = ();	# preferred to preferred
	my @freq_U_U = ();	# unpreferred to unpreferred
	my @freqS_AT_GC = (); # changes AT to GC
	my @freqS_GC_AT = (); # changes GC to AT
	my @freqS_AT_AT = (); # changes AT to AT
	my @freqS_GC_GC = (); # changes GC to GC
	my @freqR_AT_GC = (); # changes AT to GC
	my @freqR_GC_AT = (); # changes GC to AT
	my @freqR_AT_AT = (); # changes AT to AT
	my @freqR_GC_GC = (); # changes GC to GC
	my @freqS_Ts = (); # synonymous transitions
	my @freqS_Tv = (); # synonymous transversions
	
	
	$numseqs=$PR=$PS=$DR=$DS=0;

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
		# print $_,"\n";
		$numseqs++;
		}
	
	chomp $file;	
	print "\nlocus ", $file, " numseqs: ", $numseqs, "\n";
	print "outgroup used: ", $sequence_names[0], "\n";

#assign information to variable codons
	my $seqlen = length($data[0]);
#	print "sequence length before chop: $seqlen\n";
	for ($i=0; $i< ($seqlen%3); $i++)
	{
	for $k (0..(scalar(@data)-1))
	   {
	   chop ($data[$k]);
       }
	}
	my $seqlen = length($data[0]);
#	print "sequence length after chop: $seqlen\n";


# assign each codon an amino-acid and a default state NNNNNNNNN
	my @codon = ();

	$ind=$j=$k=0;

	for($j=0; $j < $seqlen; $j+=3)
		{
		for ($ind=0; $ind < $numseqs; $ind++)
			{   
			$k=(($j+3)/3)-1;
			$codon[$ind][$k] = (substr($data[$ind],$j,3));
	   	    $aa[$ind][$k]=codon2aa($codon[$ind][$k]); 
			$codon_bias[$ind][$k]=codonbias($codon[$ind][$k]); 
			if ($aa[$ind][$k] eq '_')
				{
				print "STOP codon found at position $j individual $sequence_names[$ind] in locus $file \n";
				}
			if (substr($data[$ind],$j,3) =~ tr/://)  # if no gaps returs a zero (symbol for gap "-")
			    {
			    $codon_processed[$ind][$k] = 'x';
			    }
			else 
				{
				$codon_processed[$ind][$k] = 'NNNNNNNNN';
				}
			# print "codon: ", $codon[$ind][$k], "\n";
            #  print "aa: ", $aa[$ind][$k], "\n";
			#  print "codon processed: ", $codon_processed[$ind][$k], "\n";
			# print "codon state: ", $codon_bias[$ind][$k], "\n";
			}
		 }
	
     # print "first codon, first indiv: ", $codon[0][0], "\n";
     #  print "first aa, first indiv: ", $aa[0][0], "\n";

# we have now read in the data into 4 arrays 
# codon[$ind][$k] = indiv i, codon k
# aa[$ind][$k] = corresponding amino acid for indiv i, codon k
# codon_processed[$ind][$k] = indiv i, codon k which contains default stage 'NNNNNNNNN' if there is no gap, and 'x' if there is a gap
# codon_bias[$ind][$k] = indiv i, codon k which is either preferred or non-preferred

# initialize arrays for polymorphisms ($numseqs+1 for divergence at position $numseqs+1)
for ($ind=0; $ind<($numseqs+1); $ind++)
   {
	$poly_freq_Syn[$ind]=0;
	$poly_freq_Rep[$ind]=0;
	$poly_freq_Syn_temp[$ind]=0;
	$poly_freq_Rep_temp[$ind]=0;
	$freq_P_U[$ind]=0;
	$freq_U_P[$ind]=0;
	$freq_P_P[$ind]=0;
	$freq_U_U[$ind]=0;
	$freq_P_U_temp[$ind]=0;
	$freq_U_P_temp[$ind]=0;
	$freq_P_P_temp[$ind]=0;
	$freq_U_U_temp[$ind]=0;
	$freqS_AT_GC[$ind]=0;	$freqS_GC_AT[$ind]=0;	$freqS_AT_AT[$ind]=0;	$freqS_GC_GC[$ind]=0;
	$freqR_AT_GC[$ind]=0;	$freqR_GC_AT[$ind]=0;	$freqR_AT_AT[$ind]=0;	$freqR_GC_GC[$ind]=0;
	$freqS_AT_GC_temp[$ind]=0;	$freqS_GC_AT_temp[$ind]=0;	$freqS_AT_AT_temp[$ind]=0;	$freqS_GC_GC_temp[$ind]=0;
	$freqR_AT_GC_temp[$ind]=0;	$freqR_GC_AT_temp[$ind]=0;	$freqR_AT_AT_temp[$ind]=0;	$freqR_GC_GC_temp[$ind]=0;
	$freqS_Ts[$ind]=0; $freqS_Tv[$ind]=0;
   }

$no_syn_codons=$no_rep_codons=0; $FOP=$FNOP=0; $no_syn_fourfold_codons=0; $no_syn_fourfold_div=0;
$third_pos_count=$GC_three=0; $third_pos_count_invariant=$GC_three_invariant=0;

# first task, let take polymorphic codons, and feed them to subrouting "codon_processor" if there is no gap
#**********************************************************************************************************************************************************
# big loop through all codons of a locus
for ($pos=0; $pos<($seqlen/3); $pos++)
	{
   
   # check for a gap  - only proceed if there is no gap
    $switch_gap=0;
    for ($i=0; $i<$numseqs; $i++)
       {
	   if ($aa[$i][$pos] eq 'gap')
		  {
		  $switch_gap=1;
		  }
	    }
		
  if ($switch_gap==0)		
	{
	#first count the number of synoymous and replacement sites if there is no gap
		for ($x=0; $x<$numseqs; $x++)
			{
	        $num_syn[$x]=countSyn($codon[$x][$pos]);
			$num_rep[$x]=3-$num_syn[$x];
			$no_syn_codons=$no_syn_codons+$num_syn[$x];
			$no_rep_codons=$no_rep_codons+$num_rep[$x];
			if    (codonbias($codon[$x][$pos]) eq 'P')	{$FOP++;}
			elsif (codonbias($codon[$x][$pos]) eq 'U')	{$FNOP++;}
			if (codon_fourfold($codon[$x][$pos]) == 4) {$no_syn_fourfold_codons=$no_syn_fourfold_codons+$num_syn[$x];}
			}
			
	# count GC3
		for ($x=0; $x<$numseqs; $x++)
			{
			$third_pos_count++;	
			if ((substr(($codon[$x][$pos]), 2, 1) eq 'G') or (substr(($codon[$x][$pos]), 2, 1) eq 'C')) {$GC_three++;}
			}
			
	
	# first count the number of unique codons at codon[individual][$pos] & check for codons that have divergent site 
	# that are polymorphic in the ingroup -- i.e. a base that is unique to the outgroup ($equal_bases_posX == 0) AND has
	# more than two mutations (i.e. it must be polymorphic in the ingroup - ($mut_posX > 2) at either of the three codon positions
		$x=0;
		@codon_array = ();
		@pos1_array = ();
		@pos2_array = ();
		@pos3_array = ();
		
		for ($x=0; $x<$numseqs; $x++)
			{
			$codon_array[$x]=$codon[$x][$pos];
			$pos1_array [$x]=(substr($codon[$x][$pos],0,1));
			$pos2_array [$x]=(substr($codon[$x][$pos],1,1));
			$pos3_array [$x]=(substr($codon[$x][$pos],2,1));
			# print "$codon_array[$x]\n";
			# print "$pos1_array[$x]\t$pos2_array[$x]\t$pos3_array[$x]\n";
			}
		@sorted_codon_array = sort @codon_array;
		@sorted_pos1_array = sort @pos1_array;
		@sorted_pos2_array = sort @pos2_array;
		@sorted_pos3_array = sort @pos3_array;

	#	print join("-", @sorted_pos1_array), "\n";
	#	print join("-", @sorted_pos2_array), "\n";
	#	print join("-", @sorted_pos3_array), "\n";
	
		$unique_codons=1;
		$equal_outgroup=0;
		$x=0;
		$mut_pos1=$mut_pos2=$mut_pos3=1;
		$equal_bases_pos1=$equal_bases_pos2=$equal_bases_pos3=0;
        @codons_unique =();
		push @codons_unique, $sorted_codon_array[$x];
		for ($x=0; $x<$numseqs-1; $x++)
			{
			if ($sorted_codon_array[$x] ne $sorted_codon_array[$x+1])
				{$unique_codons++;
				push @codons_unique, $sorted_codon_array[$x+1];}
			if ($sorted_pos1_array[$x] ne $sorted_pos1_array[$x+1])
				{$mut_pos1++;}
			if ($sorted_pos2_array[$x] ne $sorted_pos2_array[$x+1])	
				{$mut_pos2++;}
			if ($sorted_pos3_array[$x] ne $sorted_pos3_array[$x+1])	
				{$mut_pos3++;}

			if ($codon_array[0] eq $codon_array[$x+1])	
			    {$equal_outgroup++;}
			if ((substr($codon_array[0],0,1)) eq (substr($codon_array[$x+1],0,1)))	
			   {$equal_bases_pos1++;}
			if ((substr($codon_array[0],1,1)) eq (substr($codon_array[$x+1],1,1)))	
			   {$equal_bases_pos2++;}
			if ((substr($codon_array[0],2,1)) eq (substr($codon_array[$x+1],2,1)))	
			   {$equal_bases_pos3++;}
			}
# print "unique_codons: $unique_codons, equal to outgroup: $equal_outgroup\n";
# print " mut pos1, pos2, pos3    $mut_pos1, $mut_pos2, $mut_pos3\n";
# print " equal_bases pos1, pos2, pos3   $equal_bases_pos1, $equal_bases_pos2, $equal_bases_pos3\n"; 

# $complex = 'nein';
# if ((($equal_bases_pos1 == 0) and ($mut_pos1 > 2)) or (($equal_bases_pos2 == 0) and ($mut_pos2 > 2)) or (($equal_bases_pos3 == 0) and ($mut_pos3 > 2)))
#   {
#   $complex = 'ja';
#   }

# print "complex : $complex \n";

if ($unique_codons==1) # if monomorph -- get GC3 monomorphic
 	{
	for ($x=0; $x<$numseqs; $x++)
		{
		$third_pos_count_invariant++;	
		if ((substr(($codon[$x][$pos]), 2, 1) eq 'G') or (substr(($codon[$x][$pos]), 2, 1) eq 'C')) {$GC_three_invariant++;}
		}
	}


if ($unique_codons==2) # send data to subroutine codon_processor and get codon_processed 
 	{
	$complex = 'nein';
	@codon_array_processed = codon_processor(@codon_array, $complex);          # if codon is complicated ($complex = ja) -- return a long list	
    # print "# of entries in processed array: ", (scalar(@codon_array_processed)),"\n";

#	print "$codon_array[$pos]\t$codon_array_processed[$pos]\n";
	for ($ind=0; $ind<$numseqs; $ind++)
		{
		$codon_processed[$ind][$pos]= $codon_array_processed[$ind];
		# print "$codon_processed[$ind][$pos]\t"; 
	    }   

    # assing polymorphism and divergence
	@position1 = (); @position2 = (); @position3 = ();
	for ($ind=0; $ind<$numseqs; $ind++)
		{
        push @position1, join ("", substr($codon_array_processed[$ind], 0, 3));
		push @position2, join ("", substr($codon_array_processed[$ind], 3, 3));
		push @position3, join ("", substr($codon_array_processed[$ind], 6, 3));
		}
		
	$pos[0] = join ("", @position1);
	$pos[1] = join ("", @position2);
	$pos[2] = join ("", @position3);

#  	print "\n";
#	print "Poly at position: $pos\n";
#  	print "pos1: ", $pos[0], "\n";
#  	print "pos2: ", $pos[1], "\n";
#  	print "pos3: ", $pos[2], "\n";    
	
	@freq_NNN = (0,0,0); @freq_R = (0,0,0); @freq_S = (0,0,0);
	@freq_GR = (0,0,0); @freq_AR = (0,0,0); @freq_TR = (0,0,0); @freq_CR = (0,0,0);
	@freq_GS = (0,0,0); @freq_AS = (0,0,0); @freq_TS = (0,0,0); @freq_CS = (0,0,0);
	@freq_P = (0,0,0); @freq_U = (0,0,0); @freq_D = (0,0,0); @freq_B = (0,0,0);
	
	for $i (0..2)
		{
		while ($pos[$i] =~ /NNN/g){$freq_NNN[$i]++;}
		while ($pos[$i] =~ /S/g){$freq_S[$i]++;}
		while ($pos[$i] =~ /R/g){$freq_R[$i]++}
		while ($pos[$i] =~ /GR/g){$freq_GR[$i]++;}
		while ($pos[$i] =~ /AR/g){$freq_AR[$i]++;}
		while ($pos[$i] =~ /TR/g){$freq_TR[$i]++;}
		while ($pos[$i] =~ /CR/g){$freq_CR[$i]++;}
		while ($pos[$i] =~ /GS/g){$freq_GS[$i]++;}
		while ($pos[$i] =~ /AS/g){$freq_AS[$i]++;}
		while ($pos[$i] =~ /TS/g){$freq_TS[$i]++;}
		while ($pos[$i] =~ /CS/g){$freq_CS[$i]++;}
		while ($pos[$i] =~ /P/g){$freq_P[$i]++;}
		while ($pos[$i] =~ /U/g){$freq_U[$i]++;}
		while ($pos[$i] =~ /D/g){$freq_D[$i]++;}
		while ($pos[$i] =~ /B/g){$freq_B[$i]++;}
		}

#	 print "frequencies NNN: ", join (" ", @freq_NNN), "\t", "Syn: ", join (" ", @freq_S), "\t" , "Rep: ", join (" ", @freq_R),"\n";
#	 print "frequencies GR: ", join (" ", @freq_GR), "\t", "AR: ", join (" ", @freq_AR), "\t" , "TR: ", join (" ", @freq_TR), "\t", "CR: ", join (" ", @freq_CR), "\n";
#	 print "frequencies GS: ", join (" ", @freq_GS), "\t", "AS: ", join (" ", @freq_AS), "\t" , "TS: ", join (" ", @freq_TS), "\t", "CS: ", join (" ", @freq_CS), "\n";


# print "amino-acid: ", $aa_pos, " freq aa outgroup: $freq_aa_OG \n";
# print "codon_states: ", $codon_state_pos, "  freq P: $freq_P, freq U: $freq_U \n"; 


	for $i (0..2)	
	   {
	   $nuc_OG=substr($codon[0][$pos],$i,1);
	   if ($freq_NNN[$i] == $numseqs)   # monomorphic
	     {
		  # print "monomorphic\n";
		 }
	   elsif ($freq_NNN[$i] ==1)  # divergence
		  {
		   if ($freq_S[$i]>1)
		      {
			  $poly_freq_Syn[$numseqs]++;
			  #if (codon_fourfold($codon[0][$pos]) == 4) {$no_syn_fourfold_div++;}	#divergence at 4-fold degenerate site
			  if ($freq_P[$i]>1) {$freq_P_P[$numseqs]++;}
			  if ($freq_U[$i]>1) {$freq_U_U[$numseqs]++;}
			  if ($freq_D[$i]>1) {$freq_P_U[$numseqs]++;}
			  if ($freq_B[$i]>1) {$freq_U_P[$numseqs]++;}

			  if (($freq_GS[$i]>1) or ($freq_CS[$i]>1))  {if (($nuc_OG eq 'G') or ($nuc_OG eq 'C')) {$freqS_GC_GC[$numseqs]++;} else {$freqS_AT_GC[$numseqs]++;}}
			  if (($freq_AS[$i]>1) or ($freq_TS[$i]>1))  {if (($nuc_OG eq 'A') or ($nuc_OG eq 'T')) {$freqS_AT_AT[$numseqs]++;} else {$freqS_GC_AT[$numseqs]++;}}
			  

			  if ($nuc_OG eq 'A')   {if    ($freq_GS[$i]>1) {$freqS_Ts[$numseqs]++;} elsif ($freq_CS[$i]>1) {$freqS_Tv[$numseqs]++;}  elsif ($freq_TS[$i]>1) {$freqS_Tv[$numseqs]++;} }
			  if ($nuc_OG eq 'T')   {if    ($freq_CS[$i]>1) {$freqS_Ts[$numseqs]++;} elsif ($freq_AS[$i]>1) {$freqS_Tv[$numseqs]++;}  elsif ($freq_GS[$i]>1) {$freqS_Tv[$numseqs]++;} }
			  if ($nuc_OG eq 'G')   {if    ($freq_AS[$i]>1) {$freqS_Ts[$numseqs]++;} elsif ($freq_TS[$i]>1) {$freqS_Tv[$numseqs]++;}  elsif ($freq_CS[$i]>1) {$freqS_Tv[$numseqs]++;} }
			  if ($nuc_OG eq 'C')   {if    ($freq_TS[$i]>1) {$freqS_Ts[$numseqs]++;} elsif ($freq_AS[$i]>1) {$freqS_Tv[$numseqs]++;}  elsif ($freq_GS[$i]>1) {$freqS_Tv[$numseqs]++;} }
			  
			#   print "one synonymous divergence\n";
			  }
		   elsif ($freq_R[$i]>1)
		      {
			  $poly_freq_Rep[$numseqs]++;
			  
			  if (($freq_GR[$i]>1) or ($freq_CR[$i]>1))  {if (($nuc_OG eq 'G') or ($nuc_OG eq 'C')) {$freqR_GC_GC[$numseqs]++;} else {$freqR_AT_GC[$numseqs]++;}}
			  if (($freq_AR[$i]>1) or ($freq_TR[$i]>1))  {if (($nuc_OG eq 'A') or ($nuc_OG eq 'T')) {$freqR_AT_AT[$numseqs]++;} else {$freqR_GC_AT[$numseqs]++;}}
			  # print "one replacement divergence\n";
			  }
	       }
	  elsif (($freq_NNN[$i] >1) and ($freq_NNN[$i] < $numseqs))	# polymorphic
           {
		   
		 #    print "polymorphic\n";
			$poly_freq_Rep[$freq_GR[$i]]++;
		    $poly_freq_Rep[$freq_AR[$i]]++;
		    $poly_freq_Rep[$freq_TR[$i]]++;
		    $poly_freq_Rep[$freq_CR[$i]]++;
			$poly_freq_Syn[$freq_GS[$i]]++;
			$poly_freq_Syn[$freq_AS[$i]]++;
			$poly_freq_Syn[$freq_TS[$i]]++;
			$poly_freq_Syn[$freq_CS[$i]]++;
			$freq_P_P[$freq_P[$i]]++;
			$freq_U_U[$freq_U[$i]]++;
			$freq_P_U[$freq_D[$i]]++;
			$freq_U_P[$freq_B[$i]]++;
			
			if (($nuc_OG eq 'G') or ($nuc_OG eq 'C'))	{$freqS_GC_AT[$freq_AS[$i]]++;	$freqS_GC_AT[$freq_TS[$i]]++;	$freqS_GC_GC[$freq_GS[$i]]++;  $freqS_GC_GC[$freq_CS[$i]]++;  
														 $freqR_GC_AT[$freq_AR[$i]]++;	$freqR_GC_AT[$freq_TR[$i]]++;	$freqR_GC_GC[$freq_GR[$i]]++;  $freqR_GC_GC[$freq_CR[$i]]++;}
			if (($nuc_OG eq 'A') or ($nuc_OG eq 'T'))	{$freqS_AT_AT[$freq_AS[$i]]++;	$freqS_AT_AT[$freq_TS[$i]]++;	$freqS_AT_GC[$freq_GS[$i]]++;  $freqS_AT_GC[$freq_CS[$i]]++;  
														 $freqR_AT_AT[$freq_AR[$i]]++;	$freqR_AT_AT[$freq_TR[$i]]++;	$freqR_AT_GC[$freq_GR[$i]]++;  $freqR_AT_GC[$freq_CR[$i]]++;}

			if ($nuc_OG eq 'A')   { $freqS_Ts[$freq_GS[$i]]++; $freqS_Tv[$freq_CS[$i]]++; $freqS_Tv[$freq_TS[$i]]++; }
			if ($nuc_OG eq 'T')   { $freqS_Ts[$freq_CS[$i]]++; $freqS_Tv[$freq_AS[$i]]++; $freqS_Tv[$freq_GS[$i]]++; }
			if ($nuc_OG eq 'G')   { $freqS_Ts[$freq_AS[$i]]++; $freqS_Tv[$freq_TS[$i]]++; $freqS_Tv[$freq_CS[$i]]++; }
			if ($nuc_OG eq 'C')   { $freqS_Ts[$freq_TS[$i]]++; $freqS_Tv[$freq_AS[$i]]++; $freqS_Tv[$freq_GS[$i]]++; }
			
								 
	       }  # polymorphic
     } # through all 3 codon positions

	 
# print "polytable  Synonymous: ", join ("-", @poly_freq_Syn), "\n";
# print "polytable Replacement: ", join ("-", @poly_freq_Rep), "\n";
# print "polytable  U -> P: ", join ("-", @freq_U_P), "\n";
# print "polytable  P -> U: ", join ("-", @freq_P_U), "\n";
# print "polytable  P -> P: ", join ("-", @freq_P_P), "\n";
# print "polytable  U -> U: ", join ("-", @freq_U_U), "\n";
# print  "\n";
}   # end of if ($unique_codons==2)	

#****************************     

 if ($unique_codons>2)	
	{
	$complex = 'ja';
	@codon_array_processed = codon_processor(@codon_array, $complex);          # if codon is complicated ($complex = ja) -- return a long list	
  
	print OUT_DIFF "\nlocus ", $file, "  COMPLEX CODON at positon " ,$pos, "\n";

	
	# print "# of entries in processed array: ", (scalar(@codon_array_processed)),"\n";
	@codon_multipath_processed =();
	$num_poss_comb=(scalar(@codon_array_processed)/$numseqs);  # number of possible comparisons returned
	@NNN_codon=();
	$int_counter=0;
	for ($i=0; $i<$num_poss_comb; $i++)
	   {
	   for ($ind=0; $ind<$numseqs; $ind++)
			{
		    $codon_multipath_processed[$i][$ind]= $codon_array_processed[$int_counter];
			if ($codon_multipath_processed[$i][$ind] eq 'NNNNNNNNN')
				{
				$NNN_codon[$i]=$codon_array[$ind];
				}
		#	print "internal counter: $int_counter\n";
			$int_counter++;
			}   
          print OUT_DIFF join ("-", @{$codon_multipath_processed[$i]}), "\n";
		}
	# print "NNNNNNNNN codons: ", join ("  ", @NNN_codon), "\n";
	
	@path_min_Syn = (); @path_min_Rep = ();
	@path_min_P_U = (); @path_min_U_P = (); @path_min_P_P = (); @path_min_U_U = ();
	@path_minS_AT_GC = ();	@path_minS_AT_AT = ();	@path_minS_GC_GC = ();	@path_minS_GC_AT = ();
	@path_minR_AT_GC = ();	@path_minR_AT_AT = ();	@path_minR_GC_GC = ();	@path_minR_GC_AT = ();
	@path_minS_Ts = (); @path_minS_Tv = ();
	
	for ($ind=0; $ind<($numseqs+1); $ind++) # initialize path_min with max. no of mutations
		   {
		   $path_min_Syn[$ind]=($mut_pos3+$mut_pos2+$mut_pos1);
		   $path_min_Rep[$ind]=($mut_pos3+$mut_pos2+$mut_pos1);
		   }
	   
#	print "Syn min: $no_polyS_min, $no_divS_min \n";
#    print "Rep min: $no_polyR_min, $no_divR_min \n";

    #********** have $num_poss_comb and have to find the shortes one ************************************************************ 
    @poly_freq_Syn_multipaths = (); @poly_freq_Rep_multipaths = (); # to temoporary store all paths and find shortest
    # assing polymorphism and divergence
    for ($j=0; $j<$num_poss_comb; $j++) # loop through all possible paths
	    {
	    @poly_freq_Rep_temp=(); @poly_freq_Syn_temp=(); 
		@freq_P_U_temp = (); @freq_U_P_temp = (); @freq_P_P_temp = (); @freq_U_U_temp = ();
	    
		@position0 = (); @position1 = (); @position2 = ();

	   for ($ind=0; $ind<($numseqs+1); $ind++) # initialize temoporary arrays
		   {
		   $poly_freq_Syn_temp[$ind]=0;
		   $poly_freq_Rep_temp[$ind]=0;
		   $freq_P_U_temp[$ind]=0;
		   $freq_U_P_temp[$ind]=0;
		   $freq_P_P_temp[$ind]=0;
		   $freq_U_U_temp[$ind]=0;
		   $freqS_GC_AT_temp[$ind]=0;	$freqS_GC_GC_temp[$ind]=0;	$freqS_AT_AT_temp[$ind]=0;	$freqS_AT_GC_temp[$ind]=0;
		   $freqR_GC_AT_temp[$ind]=0;	$freqR_GC_GC_temp[$ind]=0;	$freqR_AT_AT_temp[$ind]=0;	$freqR_AT_GC_temp[$ind]=0;
		   $freqS_Ts_temp[$ind]=0;	$freqS_Tv_temp[$ind]=0;
		   }

	  for ($ind=0; $ind<$numseqs; $ind++)
		  {
          push @position0, join ("", substr($codon_multipath_processed[$j][$ind], 0, 3));
		  push @position1, join ("", substr($codon_multipath_processed[$j][$ind], 3, 3));
		  push @position2, join ("", substr($codon_multipath_processed[$j][$ind], 6, 3));
		  }

	   $sizeIG= (length $pos[0])-3;
	   $pos[0] = join ("", @position0);
	   $posOG[0] = (substr($codon_multipath_processed[$j][0], 0, 3));
	   $posIG[0] = (substr($pos[0], 3, $sizeIG));
	   $pos[1] = join ("", @position1);
	   $posOG[1] = (substr($codon_multipath_processed[$j][0], 3, 3));
	   $posIG[1] = (substr($pos[1], 3, $sizeIG));
	   $pos[2] = join ("", @position2);
       $posOG[2] = (substr($codon_multipath_processed[$j][0], 6, 3));
	   $posIG[2] = (substr($pos[2], 3, $sizeIG));
	
	  print OUT_DIFF "\nPath $j \n";
	  print OUT_DIFF "pos1: ", $pos[0],"\t outgroup ", $posOG[0],"\t ingroup ", $posIG[0], "\n";
	  print OUT_DIFF "pos2: ", $pos[1],"\t outgroup ", $posOG[1],"\t ingroup ", $posIG[1], "\n";
	  print OUT_DIFF "pos3: ", $pos[2],"\t outgroup ", $posOG[2],"\t ingroup ", $posIG[2], "\n";    
	
	  @freq_NNN = (0,0,0); @freq_R = (0,0,0); @freq_S = (0,0,0); @freq_OG = (0,0,0); @freq_MCC = (0,0,0); @freq_NNN_IG = (0,0,0);
	  @freq_GR = (0,0,0); @freq_AR = (0,0,0); @freq_TR = (0,0,0); @freq_CR = (0,0,0);
	  @freq_GR_IG = (0,0,0); @freq_AR_IG = (0,0,0); @freq_TR_IG = (0,0,0); @freq_CR_IG = (0,0,0);
	  @freq_GS = (0,0,0); @freq_AS = (0,0,0); @freq_TS = (0,0,0); @freq_CS = (0,0,0);
	  @freq_GS_IG = (0,0,0); @freq_AS_IG = (0,0,0); @freq_TS_IG = (0,0,0); @freq_CS_IG = (0,0,0);

	  @freq_P = (0,0,0); @freq_U = (0,0,0); @freq_D = (0,0,0); @freq_B = (0,0,0);
	  @freq_GSP = (0,0,0); @freq_GSU = (0,0,0); @freq_GSD = (0,0,0); @freq_GSB = (0,0,0);
	  @freq_ASP = (0,0,0); @freq_ASU = (0,0,0); @freq_ASD = (0,0,0); @freq_ASB = (0,0,0);
	  @freq_TSP = (0,0,0); @freq_TSU = (0,0,0); @freq_TSD = (0,0,0); @freq_TSB = (0,0,0);
	  @freq_CSP = (0,0,0); @freq_CSU = (0,0,0); @freq_CSD = (0,0,0); @freq_CSB = (0,0,0);
	  @most_common_codon =(0,0,0);


	  for $i (0..2)
	    {
	    while ($pos[$i] =~ /$posOG[$i]/g){$freq_OG[$i]++;}
		while ($pos[$i] =~ /S/g){$freq_S[$i]++;}
		while ($pos[$i] =~ /R/g){$freq_R[$i]++;}
		
		while ($pos[$i] =~ /NNN/g){$freq_NNN[$i]++;}
		while ($pos[$i] =~ /GR/g){$freq_GR[$i]++;}
		while ($pos[$i] =~ /AR/g){$freq_AR[$i]++;}
		while ($pos[$i] =~ /TR/g){$freq_TR[$i]++;}
		while ($pos[$i] =~ /CR/g){$freq_CR[$i]++;}
		while ($pos[$i] =~ /GS/g){$freq_GS[$i]++;}
		while ($pos[$i] =~ /AS/g){$freq_AS[$i]++;}
		while ($pos[$i] =~ /TS/g){$freq_TS[$i]++;}
		while ($pos[$i] =~ /CS/g){$freq_CS[$i]++;}

		while ($posIG[$i] =~ /NNN/g){$freq_NNN_IG[$i]++;}
		while ($posIG[$i] =~ /GR/g){$freq_GR_IG[$i]++;}
		while ($posIG[$i] =~ /AR/g){$freq_AR_IG[$i]++;}
		while ($posIG[$i] =~ /TR/g){$freq_TR_IG[$i]++;}
		while ($posIG[$i] =~ /CR/g){$freq_CR_IG[$i]++;}
		while ($posIG[$i] =~ /GS/g){$freq_GS_IG[$i]++;}
		while ($posIG[$i] =~ /AS/g){$freq_AS_IG[$i]++;}
		while ($posIG[$i] =~ /TS/g){$freq_TS_IG[$i]++;}
		while ($posIG[$i] =~ /CS/g){$freq_CS_IG[$i]++;}
		
		while ($pos[$i] =~ /P/g)	{$freq_P[$i]++;}	while ($pos[$i] =~ /U/g)	{$freq_U[$i]++;}	while ($pos[$i] =~ /D/g)	{$freq_D[$i]++;}	while ($pos[$i] =~ /B/g)	{$freq_B[$i]++;}
		while ($pos[$i] =~ /GSP/g){$freq_GSP[$i]++;}	while ($pos[$i] =~ /GSU/g){$freq_GSU[$i]++;}	while ($pos[$i] =~ /GSD/g){$freq_GSD[$i]++;}	while ($pos[$i] =~ /GSB/g){$freq_GSB[$i]++;}
		while ($pos[$i] =~ /ASP/g){$freq_ASP[$i]++;}	while ($pos[$i] =~ /ASU/g){$freq_ASU[$i]++;}	while ($pos[$i] =~ /ASD/g){$freq_ASD[$i]++;}	while ($pos[$i] =~ /ASB/g){$freq_ASB[$i]++;}
		while ($pos[$i] =~ /TSP/g){$freq_TSP[$i]++;}	while ($pos[$i] =~ /TSU/g){$freq_TSU[$i]++;}	while ($pos[$i] =~ /TSD/g){$freq_TSD[$i]++;}	while ($pos[$i] =~ /TSB/g){$freq_TSB[$i]++;}
		while ($pos[$i] =~ /CSP/g){$freq_CSP[$i]++;}	while ($pos[$i] =~ /CSU/g){$freq_CSU[$i]++;}	while ($pos[$i] =~ /CSD/g){$freq_CSD[$i]++;}	while ($pos[$i] =~ /CSB/g){$freq_CSB[$i]++;}

	    $temp=0;
		if($freq_NNN[$i]>$temp)	{$most_common='NNN'; $temp=$freq_NNN[$i];}
		if($freq_GR[$i]>$temp)	{$most_common='GR'; $temp=$freq_GR[$i];}
	    if($freq_AR[$i]>$temp)	{$most_common='AR';	$temp=$freq_AR[$i];}
		if($freq_TR[$i]>$temp)	{$most_common='TR'; $temp=$freq_TR[$i];}
		if($freq_CR[$i]>$temp)	{$most_common='CR'; $temp=$freq_CR[$i];}
		if($freq_GS[$i]>$temp)	{$most_common='GS'; $temp=$freq_GS[$i];}
		if($freq_AS[$i]>$temp)	{$most_common='AS'; $temp=$freq_AS[$i];}
		if($freq_TS[$i]>$temp)	{$most_common='TS'; $temp=$freq_TS[$i];}
		if($freq_CS[$i]>$temp)	{$most_common='CS'; $temp=$freq_CS[$i];}
		$freq_MC[$i]=$temp;	 
		$most_common_codon[$i] = $most_common;	
		}

	# print "frequencies NNN: ", join (" ", @freq_NNN), "\t", "Syn: ", join (" ", @freq_S), "\t" , "Rep: ", join (" ", @freq_R),"\t" , "OG: ", join (" ", @freq_OG),"\n";
	# print "frequencies GR   : ", join (" ", @freq_GR), "\t", "AR   : ", join (" ", @freq_AR), "\t" , "TR   : ", join (" ", @freq_TR), "\t", "CR   : ", join (" ", @freq_CR), "\n";
	# print "frequencies GR_IG: ", join (" ", @freq_GR_IG), "\t", "AR_IG: ", join (" ", @freq_AR_IG), "\t" , "TR_IG: ", join (" ", @freq_TR_IG), "\t", "CR_IG: ", join (" ", @freq_CR_IG), "\n";
	# print "frequencies GS   : ", join (" ", @freq_GS), "\t", "AS   : ", join (" ", @freq_AS), "\t" , "TS   : ", join (" ", @freq_TS), "\t", "CS   : ", join (" ", @freq_CS), "\n";
	# print "frequencies GS_IG: ", join (" ", @freq_GS_IG), "\t", "AS_IG: ", join (" ", @freq_AS_IG), "\t" , "TS_IG: ", join (" ", @freq_TS_IG), "\t", "CS_IG: ", join (" ", @freq_CS_IG), "\n";
    # print "most common codon: ", join (" ", @most_common_codon), "\t"," no. most common ",  join (" ", @freq_MC), "\n"; 	

#======================================
	$i=0;
	for $i (0..2)	
	   {
	   $nuc_NNN_codon=substr($NNN_codon[$j],$i,1);
	   if ($freq_NNN[$i] == $numseqs) 
		   {
		   #print "monomorphic \n";
		   }  #monomorphic
	   elsif ($posOG[$i] =~ /NNN/)  # outgroup is NNN
		   {
		   if (($posIG[$i] =~ /NNN/g))  # polymorphic, with ingroup having NNN as well
				{
				$poly_freq_Rep_temp[$freq_GR[$i]]++;
				$poly_freq_Rep_temp[$freq_AR[$i]]++; 
				$poly_freq_Rep_temp[$freq_TR[$i]]++;
				$poly_freq_Rep_temp[$freq_CR[$i]]++;
				$poly_freq_Syn_temp[$freq_GS[$i]]++;
				$poly_freq_Syn_temp[$freq_AS[$i]]++;
				$poly_freq_Syn_temp[$freq_TS[$i]]++;
				$poly_freq_Syn_temp[$freq_CS[$i]]++;
				
				$freq_U_U_temp[$freq_GSU[$i]]++;	$freq_U_U_temp[$freq_ASU[$i]]++;	$freq_U_U_temp[$freq_TSU[$i]]++;	$freq_U_U_temp[$freq_CSU[$i]]++; 
				$freq_P_P_temp[$freq_GSP[$i]]++;	$freq_P_P_temp[$freq_ASP[$i]]++;	$freq_P_P_temp[$freq_TSP[$i]]++;	$freq_P_P_temp[$freq_CSP[$i]]++; 
				$freq_U_P_temp[$freq_GSB[$i]]++;	$freq_U_P_temp[$freq_ASB[$i]]++;	$freq_U_P_temp[$freq_TSB[$i]]++;	$freq_U_P_temp[$freq_CSB[$i]]++; 
				$freq_P_U_temp[$freq_GSD[$i]]++;	$freq_P_U_temp[$freq_ASD[$i]]++;	$freq_P_U_temp[$freq_TSD[$i]]++;	$freq_P_U_temp[$freq_CSD[$i]]++; 
				
				if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {$freqS_GC_AT_temp[$freq_AS[$i]]++;	$freqS_GC_AT_temp[$freq_TS[$i]]++;	$freqS_GC_GC_temp[$freq_GS[$i]]++;  $freqS_GC_GC_temp[$freq_CS[$i]]++;  
																		 $freqR_GC_AT_temp[$freq_AR[$i]]++;	$freqR_GC_AT_temp[$freq_TR[$i]]++;	$freqR_GC_GC_temp[$freq_GR[$i]]++;  $freqR_GC_GC_temp[$freq_CR[$i]]++;}

				if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T'))	{$freqS_AT_AT_temp[$freq_AS[$i]]++;	$freqS_AT_AT_temp[$freq_TS[$i]]++;	$freqS_AT_GC_temp[$freq_GS[$i]]++;  $freqS_AT_GC_temp[$freq_CS[$i]]++;  
																		 $freqR_AT_AT_temp[$freq_AR[$i]]++;	$freqR_AT_AT_temp[$freq_TR[$i]]++;	$freqR_AT_GC_temp[$freq_GR[$i]]++;  $freqR_AT_GC_temp[$freq_CR[$i]]++;}
																		 
																		 
				if ($nuc_NNN_codon eq 'A')   { $freqS_Ts_temp[$freq_GS[$i]]++; $freqS_Tv_temp[$freq_CS[$i]]++; $freqS_Tv_temp[$freq_TS[$i]]++; }
				if ($nuc_NNN_codon eq 'T')   { $freqS_Ts_temp[$freq_CS[$i]]++; $freqS_Tv_temp[$freq_AS[$i]]++; $freqS_Tv_temp[$freq_GS[$i]]++; }
				if ($nuc_NNN_codon eq 'G')   { $freqS_Ts_temp[$freq_AS[$i]]++; $freqS_Tv_temp[$freq_TS[$i]]++; $freqS_Tv_temp[$freq_CS[$i]]++; }
				if ($nuc_NNN_codon eq 'C')   { $freqS_Ts_temp[$freq_TS[$i]]++; $freqS_Tv_temp[$freq_AS[$i]]++; $freqS_Tv_temp[$freq_GS[$i]]++; }
														 
																		 
				print OUT_DIFF "one (or more) poly - case: OG has NNN, and ingroup has NNN\n"
				}
			elsif (($posIG[$i] =~ /NNN/g) == 0)  # either fixed difference, or fixed difference and polymorphic 
				{
				if ($freq_MC[$i] == ($numseqs-1))  # fixed difference
				   {
				    if ($posIG[$i] =~ /R/g)
					   {
					   $poly_freq_Rep_temp[$numseqs]++;
						
						if (($freq_GR[$i]>1) or ($freq_CR[$i]>1))  {if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {$freqR_GC_GC_temp[$numseqs]++;} else {$freqR_AT_GC_temp[$numseqs]++;}}
						if (($freq_AR[$i]>1) or ($freq_TR[$i]>1))  {if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {$freqR_AT_AT_temp[$numseqs]++;} else {$freqR_GC_AT_temp[$numseqs]++;}}

						print OUT_DIFF "one replacement divergence - case: OG has NNN, no NNN in ingroup and monomorphic\n";
						}
				    if ($posIG[$i] =~ /S/g)
					   {
					    $poly_freq_Syn_temp[$numseqs]++;
					    #if (codon_fourfold($codon[0][$pos]) == 4) {$no_syn_fourfold_div++;}	#divergence at 4-fold degenerate site
						if ($freq_U[$i] == ($numseqs-1)) {$freq_U_U_temp[$numseqs]++;}
						if ($freq_P[$i] == ($numseqs-1)) {$freq_P_P_temp[$numseqs]++;}
						if ($freq_B[$i] == ($numseqs-1)) {$freq_U_P_temp[$numseqs]++;}
						if ($freq_D[$i] == ($numseqs-1)) {$freq_P_U_temp[$numseqs]++;}
						
						if (($freq_GS[$i]>1) or ($freq_CS[$i]>1))  {if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {$freqS_GC_GC_temp[$numseqs]++;} else {$freqS_AT_GC_temp[$numseqs]++;}}
						if (($freq_AS[$i]>1) or ($freq_TS[$i]>1))  {if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {$freqS_AT_AT_temp[$numseqs]++;} else {$freqS_GC_AT_temp[$numseqs]++;}}
						
						if ($nuc_NNN_codon eq 'A')   {if    ($freq_GS[$i]>1) {$freqS_Ts_temp[$numseqs]++;} elsif ($freq_CS[$i]>1) {$freqS_Tv_temp[$numseqs]++;}  elsif ($freq_TS[$i]>1) {$freqS_Tv_temp[$numseqs]++;} }
						if ($nuc_NNN_codon eq 'T')   {if    ($freq_CS[$i]>1) {$freqS_Ts_temp[$numseqs]++;} elsif ($freq_AS[$i]>1) {$freqS_Tv_temp[$numseqs]++;}  elsif ($freq_GS[$i]>1) {$freqS_Tv_temp[$numseqs]++;} }
						if ($nuc_NNN_codon eq 'G')   {if    ($freq_AS[$i]>1) {$freqS_Ts_temp[$numseqs]++;} elsif ($freq_TS[$i]>1) {$freqS_Tv_temp[$numseqs]++;}  elsif ($freq_CS[$i]>1) {$freqS_Tv_temp[$numseqs]++;} }
						if ($nuc_NNN_codon eq 'C')   {if    ($freq_TS[$i]>1) {$freqS_Ts_temp[$numseqs]++;} elsif ($freq_AS[$i]>1) {$freqS_Tv_temp[$numseqs]++;}  elsif ($freq_GS[$i]>1) {$freqS_Tv_temp[$numseqs]++;} }

						
						print OUT_DIFF "one synonymnous divergence - case: OG has NNN, no NNN in ingroup and monomorphic\n";
						}
					}
				else # no NNN in ingroup (i.e. divergence) and polymorphic - take most common codon and make it ancestor
					{
					if (($most_common_codon[$i] =~ /R/g))
						{
						$poly_freq_Rep_temp[$numseqs]++;
						
						if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {$freqR_GC_GC_temp[$numseqs]++;} else {$freqR_AT_GC_temp[$numseqs]++;}}
						if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {$freqR_AT_AT_temp[$numseqs]++;} else {$freqR_GC_AT_temp[$numseqs]++;}}
						print OUT_DIFF "one replacement divergence- case: OG has NNN, no NNN in ingroup (i.e. divergence) and polymorphic\n";
						}
					if (($most_common_codon[$i] =~ /S/g))
						{
						$poly_freq_Syn_temp[$numseqs]++;
						#if (codon_fourfold($most_common_codon) == 4) {$no_syn_fourfold_div++;}	#divergence at 4-fold degenerate site
						if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {$freqS_GC_GC_temp[$numseqs]++;} else {$freqS_AT_GC_temp[$numseqs]++;}}
						if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {$freqS_AT_AT_temp[$numseqs]++;} else {$freqS_GC_AT_temp[$numseqs]++;}}
						if ('GS' eq $most_common_codon[$i]) {	if($freq_GSU[$i] != 0) {$freq_U_U_temp[$numseqs]++;}
																if($freq_GSP[$i] != 0) {$freq_P_P_temp[$numseqs]++;}	
																if($freq_GSB[$i] != 0) {$freq_U_P_temp[$numseqs]++;}	
																if($freq_GSD[$i] != 0) {$freq_P_U_temp[$numseqs]++;} }
						if ('AS' eq $most_common_codon[$i]) {	if($freq_ASU[$i] != 0) {$freq_U_U_temp[$numseqs]++;}
																if($freq_ASP[$i] != 0) {$freq_P_P_temp[$numseqs]++;}	
																if($freq_ASB[$i] != 0) {$freq_U_P_temp[$numseqs]++;}	
																if($freq_ASD[$i] != 0) {$freq_P_U_temp[$numseqs]++;} }
						if ('TS' eq $most_common_codon[$i]) {	if($freq_TSU[$i] != 0) {$freq_U_U_temp[$numseqs]++;}
																if($freq_TSP[$i] != 0) {$freq_P_P_temp[$numseqs]++;}	
																if($freq_TSB[$i] != 0) {$freq_U_P_temp[$numseqs]++;}	
																if($freq_TSD[$i] != 0) {$freq_P_U_temp[$numseqs]++;} }
						if ('CS' eq $most_common_codon[$i]) {	if($freq_CSU[$i] != 0) {$freq_U_U_temp[$numseqs]++;}
																if($freq_CSP[$i] != 0) {$freq_P_P_temp[$numseqs]++;}	
																if($freq_CSB[$i] != 0) {$freq_U_P_temp[$numseqs]++;}	
																if($freq_CSD[$i] != 0) {$freq_P_U_temp[$numseqs]++;} }
																
						if ($most_common_codon[$i] =~ /A/g)   {if    ($nuc_NNN_codon eq 'G') {$freqS_Ts_temp[$numseqs]++;} elsif ($nuc_NNN_codon eq 'C') {$freqS_Tv_temp[$numseqs]++;}  elsif ($nuc_NNN_codon eq 'T') {$freqS_Tv_temp[$numseqs]++;} }
						if ($most_common_codon[$i] =~ /T/g)   {if    ($nuc_NNN_codon eq 'C') {$freqS_Ts_temp[$numseqs]++;} elsif ($nuc_NNN_codon eq 'A') {$freqS_Tv_temp[$numseqs]++;}  elsif ($nuc_NNN_codon eq 'G') {$freqS_Tv_temp[$numseqs]++;} }
						if ($most_common_codon[$i] =~ /G/g)   {if    ($nuc_NNN_codon eq 'A') {$freqS_Ts_temp[$numseqs]++;} elsif ($nuc_NNN_codon eq 'T') {$freqS_Tv_temp[$numseqs]++;}  elsif ($nuc_NNN_codon eq 'C') {$freqS_Tv_temp[$numseqs]++;} }
						if ($most_common_codon[$i] =~ /C/g)   {if    ($nuc_NNN_codon eq 'T') {$freqS_Ts_temp[$numseqs]++;} elsif ($nuc_NNN_codon eq 'A') {$freqS_Tv_temp[$numseqs]++;}  elsif ($nuc_NNN_codon eq 'G') {$freqS_Tv_temp[$numseqs]++;} }
										
																
						print OUT_DIFF "one synonymnous divergence- case: OG has NNN, no NNN in ingroup (i.e. divergence) and polymorphic\n";
						}
					if ('GR' ne $most_common_codon[$i])	{$poly_freq_Rep_temp[$freq_GR[$i]]++; print OUT_DIFF "GR freq: ", $freq_GR[$i],"\n";}
					if ('AR' ne $most_common_codon[$i])	{$poly_freq_Rep_temp[$freq_AR[$i]]++; print OUT_DIFF "AR freq: ", $freq_AR[$i],"\n";}	
					if ('TR' ne $most_common_codon[$i])	{$poly_freq_Rep_temp[$freq_TR[$i]]++; print OUT_DIFF "TR freq: ", $freq_TR[$i],"\n";}
					if ('CR' ne $most_common_codon[$i])	{$poly_freq_Rep_temp[$freq_CR[$i]]++; print OUT_DIFF "CR freq: ", $freq_CR[$i],"\n";}
					if ('GS' ne $most_common_codon[$i])	{$poly_freq_Syn_temp[$freq_GS[$i]]++;	$freq_U_U_temp[$freq_GSU[$i]]++;	$freq_P_P_temp[$freq_GSP[$i]]++;	$freq_U_P_temp[$freq_GSB[$i]]++;	$freq_P_U_temp[$freq_GSD[$i]]++;	print OUT_DIFF "GS freq: ", $freq_GS[$i],"\n";}
					if ('AS' ne $most_common_codon[$i])	{$poly_freq_Syn_temp[$freq_AS[$i]]++;	$freq_U_U_temp[$freq_ASU[$i]]++;	$freq_P_P_temp[$freq_ASP[$i]]++;	$freq_U_P_temp[$freq_ASB[$i]]++;	$freq_P_U_temp[$freq_ASD[$i]]++;	print OUT_DIFF "AS freq: ", $freq_AS[$i],"\n";}
					if ('TS' ne $most_common_codon[$i])	{$poly_freq_Syn_temp[$freq_TS[$i]]++;	$freq_U_U_temp[$freq_TSU[$i]]++;	$freq_P_P_temp[$freq_TSP[$i]]++;	$freq_U_P_temp[$freq_TSB[$i]]++;	$freq_P_U_temp[$freq_TSD[$i]]++;	print OUT_DIFF "TS freq: ", $freq_TS[$i],"\n";}
					if ('CS' ne $most_common_codon[$i])	{$poly_freq_Syn_temp[$freq_CS[$i]]++;	$freq_U_U_temp[$freq_CSU[$i]]++;	$freq_P_P_temp[$freq_CSP[$i]]++;	$freq_U_P_temp[$freq_CSB[$i]]++;	$freq_P_U_temp[$freq_CSD[$i]]++;	print OUT_DIFF "CS freq: ", $freq_CS[$i],"\n";}

					if ('GR' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqR_GC_GC_temp[$freq_GR[$i]]++;} else {$freqR_AT_GC_temp[$freq_GR[$i]]++;} }
					if ('AR' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqR_AT_AT_temp[$freq_AR[$i]]++;} else {$freqR_GC_AT_temp[$freq_AR[$i]]++;} }	
					if ('TR' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqR_AT_AT_temp[$freq_TR[$i]]++;} else {$freqR_GC_AT_temp[$freq_TR[$i]]++;} }
					if ('CR' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqR_GC_GC_temp[$freq_CR[$i]]++;} else {$freqR_AT_GC_temp[$freq_CR[$i]]++;} }
					if ('GS' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqS_GC_GC_temp[$freq_GS[$i]]++;} else {$freqS_AT_GC_temp[$freq_GS[$i]]++;} }
					if ('AS' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqS_AT_AT_temp[$freq_AS[$i]]++;} else {$freqS_GC_AT_temp[$freq_AS[$i]]++;} }
					if ('TS' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqS_AT_AT_temp[$freq_TS[$i]]++;} else {$freqS_GC_AT_temp[$freq_TS[$i]]++;} }
					if ('CS' ne $most_common_codon[$i])	{if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqS_GC_GC_temp[$freq_CS[$i]]++;} else {$freqS_AT_GC_temp[$freq_CS[$i]]++;} }

					
					if ('GS' ne $most_common_codon[$i])	{if ($most_common_codon[$i] =~ /AS/g) { $freqS_Ts_temp[$freq_GS[$i]]++;} else {$freqS_Tv_temp[$freq_GS[$i]]++;} }
					if ('AS' ne $most_common_codon[$i])	{if ($most_common_codon[$i] =~ /GS/g) { $freqS_Ts_temp[$freq_AS[$i]]++;} else {$freqS_Tv_temp[$freq_AS[$i]]++;} }
					if ('TS' ne $most_common_codon[$i])	{if ($most_common_codon[$i] =~ /CS/g) { $freqS_Ts_temp[$freq_TS[$i]]++;} else {$freqS_Tv_temp[$freq_TS[$i]]++;} }
					if ('CS' ne $most_common_codon[$i])	{if ($most_common_codon[$i] =~ /TS/g) { $freqS_Ts_temp[$freq_CS[$i]]++;} else {$freqS_Tv_temp[$freq_CS[$i]]++;} }

					}
				}
			}	
		elsif ($posOG[$i] !~ /NNN/g)  # outgroup has no NNN
			{
			   if ($freq_OG[$i]>1) # outgroup is present in sample - no divergence, increment R or S poly (depending on OG) for NNN-class
					{
					if (($posOG[$i] =~ /R/g))
						 {
						 $poly_freq_Rep_temp[$freq_NNN[$i]]++;
						 if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqR_GC_GC_temp[$freq_NNN[$i]]++;} else {$freqR_AT_GC_temp[$freq_NNN[$i]]++;}}
						 if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqR_AT_AT_temp[$freq_NNN[$i]]++;} else {$freqR_GC_AT_temp[$freq_NNN[$i]]++;}}
						 print OUT_DIFF "one replacement poly- case: OG no NNN, but present in sample\n";
						 }
					if (($posOG[$i]  =~ /S/g))
						 {
						 $poly_freq_Syn_temp[$freq_NNN[$i]]++;
						 if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqS_GC_GC_temp[$freq_NNN[$i]]++;} else {$freqS_AT_GC_temp[$freq_NNN[$i]]++;}}
						 if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqS_AT_AT_temp[$freq_NNN[$i]]++;} else {$freqS_GC_AT_temp[$freq_NNN[$i]]++;}}
						 if (($posOG[$i]  =~ /U/g))	{$freq_U_U_temp[$freq_NNN[$i]]++;}
						 if (($posOG[$i]  =~ /P/g))	{$freq_P_P_temp[$freq_NNN[$i]]++;}
						 if (($posOG[$i]  =~ /B/g))	{$freq_U_P_temp[$freq_NNN[$i]]++;}
						 if (($posOG[$i]  =~ /D/g))	{$freq_P_U_temp[$freq_NNN[$i]]++;}
					
						if ($nuc_NNN_codon eq 'A')   {if    ($posOG[$i] =~ /G/g) {$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}   }
						if ($nuc_NNN_codon eq 'T')   {if    ($posOG[$i] =~ /C/g) {$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}   }
						if ($nuc_NNN_codon eq 'G')   {if    ($posOG[$i] =~ /A/g) {$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}   }
						if ($nuc_NNN_codon eq 'C')   {if    ($posOG[$i] =~ /T/g) {$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}   }

						 print OUT_DIFF "one synonymnous poly- case: OG no NNN, but present in sample\n";
						}

					if ('GR' ne $posOG[$i])	{$poly_freq_Rep_temp[$freq_GR[$i]]++; print OUT_DIFF "GR freq: ", $freq_GR[$i],"\n";}
					if ('AR' ne $posOG[$i])	{$poly_freq_Rep_temp[$freq_AR[$i]]++; print OUT_DIFF "AR freq: ", $freq_AR[$i],"\n";}	
					if ('TR' ne $posOG[$i])	{$poly_freq_Rep_temp[$freq_TR[$i]]++; print OUT_DIFF "TR freq: ", $freq_TR[$i],"\n";}
					if ('CR' ne $posOG[$i])	{$poly_freq_Rep_temp[$freq_CR[$i]]++; print OUT_DIFF "CR freq: ", $freq_CR[$i],"\n";}
					if ('GS' ne $posOG[$i])	{$poly_freq_Syn_temp[$freq_GS[$i]]++;	$freq_U_U_temp[$freq_GSU[$i]]++;	$freq_P_P_temp[$freq_GSP[$i]]++;	$freq_U_P_temp[$freq_GSB[$i]]++;	$freq_P_U_temp[$freq_GSD[$i]]++;	print OUT_DIFF "GS freq: ", $freq_GS[$i],"\n";}
					if ('AS' ne $posOG[$i])	{$poly_freq_Syn_temp[$freq_AS[$i]]++;	$freq_U_U_temp[$freq_ASU[$i]]++;	$freq_P_P_temp[$freq_ASP[$i]]++;	$freq_U_P_temp[$freq_ASB[$i]]++;	$freq_P_U_temp[$freq_ASD[$i]]++;	print OUT_DIFF "AS freq: ", $freq_AS[$i],"\n";}
					if ('TS' ne $posOG[$i])	{$poly_freq_Syn_temp[$freq_TS[$i]]++;	$freq_U_U_temp[$freq_TSU[$i]]++;	$freq_P_P_temp[$freq_TSP[$i]]++;	$freq_U_P_temp[$freq_TSB[$i]]++;	$freq_P_U_temp[$freq_TSD[$i]]++;	print OUT_DIFF "TS freq: ", $freq_TS[$i],"\n";}
					if ('CS' ne $posOG[$i])	{$poly_freq_Syn_temp[$freq_CS[$i]]++;	$freq_U_U_temp[$freq_CSU[$i]]++;	$freq_P_P_temp[$freq_CSP[$i]]++;	$freq_U_P_temp[$freq_CSB[$i]]++;	$freq_P_U_temp[$freq_CSD[$i]]++;	print OUT_DIFF "CS freq: ", $freq_CS[$i],"\n";}

					if ('GR' ne $posOG[$i]) {if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqR_GC_GC_temp[$freq_GR[$i]]++;} else {$freqR_AT_GC_temp[$freq_GR[$i]]++;} }
					if ('AR' ne $posOG[$i])	{if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqR_AT_AT_temp[$freq_AR[$i]]++;} else {$freqR_GC_AT_temp[$freq_AR[$i]]++;} }
					if ('TR' ne $posOG[$i])	{if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqR_AT_AT_temp[$freq_TR[$i]]++;} else {$freqR_GC_AT_temp[$freq_TR[$i]]++;} }
					if ('CR' ne $posOG[$i])	{if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqR_GC_GC_temp[$freq_CR[$i]]++;} else {$freqR_AT_GC_temp[$freq_CR[$i]]++;} }
					if ('GS' ne $posOG[$i])	{if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqS_GC_GC_temp[$freq_GS[$i]]++;} else {$freqS_AT_GC_temp[$freq_GS[$i]]++;} }
					if ('AS' ne $posOG[$i])	{if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqS_AT_AT_temp[$freq_AS[$i]]++;} else {$freqS_GC_AT_temp[$freq_AS[$i]]++;} }
					if ('TS' ne $posOG[$i])	{if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqS_AT_AT_temp[$freq_TS[$i]]++;} else {$freqS_GC_AT_temp[$freq_TS[$i]]++;} }
					if ('CS' ne $posOG[$i])	{if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqS_GC_GC_temp[$freq_CS[$i]]++;} else {$freqS_AT_GC_temp[$freq_CS[$i]]++;} }
					
					if ('GS' ne $posOG[$i])	{if ($posOG[$i] =~ /A/g) { $freqS_Ts_temp[$freq_GS[$i]]++;} else {$freqS_Tv_temp[$freq_GS[$i]]++;} }
					if ('AS' ne $posOG[$i])	{if ($posOG[$i] =~ /G/g) { $freqS_Ts_temp[$freq_AS[$i]]++;} else {$freqS_Tv_temp[$freq_AS[$i]]++;} }
					if ('TS' ne $posOG[$i])	{if ($posOG[$i] =~ /C/g) { $freqS_Ts_temp[$freq_TS[$i]]++;} else {$freqS_Tv_temp[$freq_TS[$i]]++;} }
					if ('CS' ne $posOG[$i])	{if ($posOG[$i] =~ /T/g) { $freqS_Ts_temp[$freq_CS[$i]]++;} else {$freqS_Tv_temp[$freq_CS[$i]]++;} }
					
					}
			   elsif ($freq_OG[$i]==1)
					{
					if ($freq_NNN[$i] == ($numseqs-1))  # fixed difference
						{
						if ($posOG[$i] =~ /R/g)
							{
							$poly_freq_Rep_temp[$numseqs]++;
							if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqR_GC_GC_temp[$numseqs]++;} else {$freqR_AT_GC_temp[$numseqs]++;}}
							if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqR_AT_AT_temp[$numseqs]++;} else {$freqR_GC_AT_temp[$numseqs]++;}}
							print OUT_DIFF "one replacement divergence - case: OG no NNN, but sample monomorphic\n";
							}
						if ($posOG[$i] =~ /S/g)
							{
							$poly_freq_Syn_temp[$numseqs]++;
							#if (codon_fourfold($codon[0][$pos]) == 4) {$no_syn_fourfold_div++;}	#divergence at 4-fold degenerate site
							if (($posOG[$i]  =~ /U/g))	{$freq_U_U_temp[$numseqs]++;}
							if (($posOG[$i]  =~ /P/g))	{$freq_P_P_temp[$numseqs]++;}
							if (($posOG[$i]  =~ /B/g))	{$freq_U_P_temp[$numseqs]++;}
							if (($posOG[$i]  =~ /D/g))	{$freq_P_U_temp[$numseqs]++;}
							if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqS_GC_GC_temp[$numseqs]++;} else {$freqS_AT_GC_temp[$numseqs]++;}}
							if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqS_AT_AT_temp[$numseqs]++;} else {$freqS_GC_AT_temp[$numseqs]++;}}
							
							if ($nuc_NNN_codon eq 'A')  {if ($posOG[$i] =~ /G/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;} }
							if ($nuc_NNN_codon eq 'G')  {if ($posOG[$i] =~ /A/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;} }
							if ($nuc_NNN_codon eq 'C')  {if ($posOG[$i] =~ /T/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;} }
							if ($nuc_NNN_codon eq 'T')  {if ($posOG[$i] =~ /C/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;} }
							
							print OUT_DIFF "one synonymnous divergence - case: OG no NNN, but sample monomorphic\n";
							}
						}
					else # OG not in ingroup (i.e. divergence) and polymorphic - take divergence from outgroup and make most common codon ancestor
						{
						if ($posOG[$i] =~ /R/g)
							{
							$poly_freq_Rep_temp[$numseqs]++;
							if ($most_common_codon[$i] eq 'NNN') { if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C'))	{if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqR_GC_GC_temp[$numseqs]++;} else {$freqR_GC_AT_temp[$numseqs]++;}}
																	else													{if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqR_AT_AT_temp[$numseqs]++;} else {$freqR_AT_GC_temp[$numseqs]++;}}}
							else  {	if (($most_common_codon[$i]=~ /G/g) or ($most_common_codon[$i]=~ /C/g)) {if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqR_GC_GC_temp[$numseqs]++;} else {$freqR_GC_AT_temp[$numseqs]++;}}
									if (($most_common_codon[$i]=~ /A/g) or ($most_common_codon[$i]=~ /T/g)) {if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqR_AT_AT_temp[$numseqs]++;} else {$freqR_AT_GC_temp[$numseqs]++;}}}
							print OUT_DIFF "one replacement divergence - case: OG no NNN, and sample polymorphic\n";
							}
						if ($posOG[$i] =~ /S/g)
							{
							$poly_freq_Syn_temp[$numseqs]++;
							#if (codon_fourfold($most_common_codon) == 4) {$no_syn_fourfold_div++;}	#divergence at 4-fold degenerate site

							if ($most_common_codon[$i] eq 'NNN') { if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C'))	{if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqS_GC_GC_temp[$numseqs]++;} else {$freqS_GC_AT_temp[$numseqs]++;}}
																	else													{if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqS_AT_AT_temp[$numseqs]++;} else {$freqS_AT_GC_temp[$numseqs]++;}}}
							else  {	if (($most_common_codon[$i]=~ /G/g) or ($most_common_codon[$i]=~ /C/g)) {if (($posOG[$i] =~ /G/g) or ($posOG[$i] =~ /C/g)) {$freqS_GC_GC_temp[$numseqs]++;} else {$freqS_GC_AT_temp[$numseqs]++;}}
									if (($most_common_codon[$i]=~ /A/g) or ($most_common_codon[$i]=~ /T/g)) {if (($posOG[$i] =~ /A/g) or ($posOG[$i] =~ /T/g)) {$freqS_AT_AT_temp[$numseqs]++;} else {$freqS_AT_GC_temp[$numseqs]++;}}}
							if (($posOG[$i]  =~ /U/g))	{$freq_U_U_temp[$numseqs]++;}
							if (($posOG[$i]  =~ /P/g))	{$freq_P_P_temp[$numseqs]++;}
							if (($posOG[$i]  =~ /B/g))	{$freq_U_P_temp[$numseqs]++;}
							if (($posOG[$i]  =~ /D/g))	{$freq_P_U_temp[$numseqs]++;}
							
							if ($most_common_codon[$i] eq 'NNN') {	if ($nuc_NNN_codon eq 'A')  {if ($posOG[$i] =~ /G/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}
																	if ($nuc_NNN_codon eq 'G')  {if ($posOG[$i] =~ /A/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}
																	if ($nuc_NNN_codon eq 'C')  {if ($posOG[$i] =~ /T/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}
																	if ($nuc_NNN_codon eq 'T')  {if ($posOG[$i] =~ /C/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}	}
							else  {	if ($most_common_codon[$i]=~ /A/g) {if ($posOG[$i] =~ /G/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}
									if ($most_common_codon[$i]=~ /G/g) {if ($posOG[$i] =~ /A/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}
									if ($most_common_codon[$i]=~ /C/g) {if ($posOG[$i] =~ /T/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}
									if ($most_common_codon[$i]=~ /T/g) {if ($posOG[$i] =~ /C/g)  {$freqS_Ts_temp[$numseqs]++;} else {$freqS_Tv_temp[$numseqs]++;}}}

							print OUT_DIFF "one synonymnous divergence- case: OG no NNN, and sample polymorphic\n";
							}
							
						if (('GR' ne $most_common_codon[$i]) and ('GR' ne $posOG[$i])) {$poly_freq_Rep_temp[$freq_GR[$i]]++; }
						if (('AR' ne $most_common_codon[$i]) and ('AR' ne $posOG[$i])) {$poly_freq_Rep_temp[$freq_AR[$i]]++; }	
						if (('TR' ne $most_common_codon[$i]) and ('TR' ne $posOG[$i])) {$poly_freq_Rep_temp[$freq_TR[$i]]++; }
						if (('CR' ne $most_common_codon[$i]) and ('CR' ne $posOG[$i])) {$poly_freq_Rep_temp[$freq_CR[$i]]++; }
						if (('GS' ne $most_common_codon[$i]) and ('GS' ne $posOG[$i])) {$poly_freq_Syn_temp[$freq_GS[$i]]++;	$freq_U_U_temp[$freq_GSU[$i]]++;	$freq_P_P_temp[$freq_GSP[$i]]++;	$freq_U_P_temp[$freq_GSB[$i]]++;	$freq_P_U_temp[$freq_GSD[$i]]++;}
						if (('AS' ne $most_common_codon[$i]) and ('AS' ne $posOG[$i])) {$poly_freq_Syn_temp[$freq_AS[$i]]++;	$freq_U_U_temp[$freq_ASU[$i]]++;	$freq_P_P_temp[$freq_ASP[$i]]++;	$freq_U_P_temp[$freq_ASB[$i]]++;	$freq_P_U_temp[$freq_ASD[$i]]++;}
						if (('TS' ne $most_common_codon[$i]) and ('TS' ne $posOG[$i])) {$poly_freq_Syn_temp[$freq_TS[$i]]++;	$freq_U_U_temp[$freq_TSU[$i]]++;	$freq_P_P_temp[$freq_TSP[$i]]++;	$freq_U_P_temp[$freq_TSB[$i]]++;	$freq_P_U_temp[$freq_TSD[$i]]++;}
						if (('CS' ne $most_common_codon[$i]) and ('CS' ne $posOG[$i])) {$poly_freq_Syn_temp[$freq_CS[$i]]++;	$freq_U_U_temp[$freq_CSU[$i]]++;	$freq_P_P_temp[$freq_CSP[$i]]++;	$freq_U_P_temp[$freq_CSB[$i]]++;	$freq_P_U_temp[$freq_CSD[$i]]++;}
						if ('NNN' ne $most_common_codon[$i])
							{
							if ($most_common_codon[$i]=~ /R/g)	{$poly_freq_Rep_temp[$freq_NNN[$i]]++; print OUT_DIFF "NNN-R freq: ", $freq_NNN[$i],"\n";}
							if ($most_common_codon[$i]=~ /S/g)	{$poly_freq_Syn_temp[$freq_NNN[$i]]++; print OUT_DIFF "NNN-S freq: ", $freq_NNN[$i],"\n";
																	if ('GS' eq $most_common_codon[$i]) {	if($freq_GSU[$i] != 0) {$freq_U_U_temp[$freq_NNN[$i]]++;}
																											if($freq_GSP[$i] != 0) {$freq_P_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_GSB[$i] != 0) {$freq_U_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_GSD[$i] != 0) {$freq_P_U_temp[$freq_NNN[$i]]++;} }
																	if ('AS' eq $most_common_codon[$i]) {	if($freq_ASU[$i] != 0) {$freq_U_U_temp[$freq_NNN[$i]]++;}
																											if($freq_ASP[$i] != 0) {$freq_P_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_ASB[$i] != 0) {$freq_U_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_ASD[$i] != 0) {$freq_P_U_temp[$freq_NNN[$i]]++;} }
																	if ('TS' eq $most_common_codon[$i]) {	if($freq_TSU[$i] != 0) {$freq_U_U_temp[$freq_NNN[$i]]++;}
																											if($freq_TSP[$i] != 0) {$freq_P_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_TSB[$i] != 0) {$freq_U_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_TSD[$i] != 0) {$freq_P_U_temp[$freq_NNN[$i]]++;} }
																	if ('CS' eq $most_common_codon[$i]) {	if($freq_CSU[$i] != 0) {$freq_U_U_temp[$freq_NNN[$i]]++;}
																											if($freq_CSP[$i] != 0) {$freq_P_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_CSB[$i] != 0) {$freq_U_P_temp[$freq_NNN[$i]]++;}	
																											if($freq_CSD[$i] != 0) {$freq_P_U_temp[$freq_NNN[$i]]++;} }
																}
							}
							
							
						if (('GR' ne $most_common_codon[$i]) and ('GR' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C'))					{$freqR_GC_GC_temp[$freq_GR[$i]]++;} else {$freqR_AT_GC_temp[$freq_GR[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqR_GC_GC_temp[$freq_GR[$i]]++;} else {$freqR_AT_GC_temp[$freq_GR[$i]]++;}}}
						if (('AR' ne $most_common_codon[$i]) and ('AR' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T'))					{$freqR_AT_AT_temp[$freq_AR[$i]]++;} else {$freqR_GC_AT_temp[$freq_AR[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqR_AT_AT_temp[$freq_AR[$i]]++;} else {$freqR_GC_AT_temp[$freq_AR[$i]]++;}}}	
						if (('TR' ne $most_common_codon[$i]) and ('TR' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T'))					{$freqR_AT_AT_temp[$freq_TR[$i]]++;} else {$freqR_GC_AT_temp[$freq_TR[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqR_AT_AT_temp[$freq_TR[$i]]++;} else {$freqR_GC_AT_temp[$freq_TR[$i]]++;}}}
						if (('CR' ne $most_common_codon[$i]) and ('CR' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C'))					{$freqR_GC_GC_temp[$freq_CR[$i]]++;} else {$freqR_AT_GC_temp[$freq_CR[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqR_GC_GC_temp[$freq_CR[$i]]++;} else {$freqR_AT_GC_temp[$freq_CR[$i]]++;}}}
						if (('GS' ne $most_common_codon[$i]) and ('GS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C'))					{$freqS_GC_GC_temp[$freq_GS[$i]]++;} else {$freqS_AT_GC_temp[$freq_GS[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqS_GC_GC_temp[$freq_GS[$i]]++;} else {$freqS_AT_GC_temp[$freq_GS[$i]]++;}}}
						if (('AS' ne $most_common_codon[$i]) and ('AS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T'))					{$freqS_AT_AT_temp[$freq_AS[$i]]++;} else {$freqS_GC_AT_temp[$freq_AS[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqS_AT_AT_temp[$freq_AS[$i]]++;} else {$freqS_GC_AT_temp[$freq_AS[$i]]++;}}}
						if (('TS' ne $most_common_codon[$i]) and ('TS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T'))					{$freqS_AT_AT_temp[$freq_TS[$i]]++;} else {$freqS_GC_AT_temp[$freq_TS[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /A/g) or ($most_common_codon[$i] =~ /T/g)) {$freqS_AT_AT_temp[$freq_TS[$i]]++;} else {$freqS_GC_AT_temp[$freq_TS[$i]]++;}}}
						if (('CS' ne $most_common_codon[$i]) and ('CS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C'))					{$freqS_GC_GC_temp[$freq_CS[$i]]++;} else {$freqS_AT_GC_temp[$freq_CS[$i]]++;}}
																						else									{ if (($most_common_codon[$i] =~ /G/g) or ($most_common_codon[$i] =~ /C/g)) {$freqS_GC_GC_temp[$freq_CS[$i]]++;} else {$freqS_AT_GC_temp[$freq_CS[$i]]++;}}}
						if ('NNN' ne $most_common_codon[$i])	{if ($most_common_codon[$i]=~ /R/g)	{if (($most_common_codon[$i]=~ /G/g) or ($most_common_codon[$i]=~ /C/g)) { if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {$freqR_GC_GC_temp[$freq_NNN[$i]]++;} else {$freqR_GC_AT_temp[$freq_NNN[$i]]++;}}
																									 else																	 { if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {$freqR_AT_AT_temp[$freq_NNN[$i]]++;} else {$freqR_AT_GC_temp[$freq_NNN[$i]]++;}}}
																 if ($most_common_codon[$i]=~ /S/g)	{if (($most_common_codon[$i]=~ /G/g) or ($most_common_codon[$i]=~ /C/g)) { if (($nuc_NNN_codon eq 'G') or ($nuc_NNN_codon eq 'C')) {$freqS_GC_GC_temp[$freq_NNN[$i]]++;} else {$freqS_GC_AT_temp[$freq_NNN[$i]]++;}}
																									 else																	 { if (($nuc_NNN_codon eq 'A') or ($nuc_NNN_codon eq 'T')) {$freqS_AT_AT_temp[$freq_NNN[$i]]++;} else {$freqS_AT_GC_temp[$freq_NNN[$i]]++;}}}
																}

# i am here 
						if (('GS' ne $most_common_codon[$i]) and ('GS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if ($nuc_NNN_codon eq 'A')			{$freqS_Ts_temp[$freq_GS[$i]]++;}  else {$freqS_Tv_temp[$freq_GS[$i]]++;}}
																						else									{ if ($most_common_codon[$i] =~ /A/g)	{$freqS_Ts_temp[$freq_GS[$i]]++;}  else {$freqS_Tv_temp[$freq_GS[$i]]++;}}}
						if (('AS' ne $most_common_codon[$i]) and ('AS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if ($nuc_NNN_codon eq 'G') 			{$freqS_Ts_temp[$freq_AS[$i]]++;}  else {$freqS_Tv_temp[$freq_AS[$i]]++;}}
																						else									{ if ($most_common_codon[$i] =~ /G/g)	{$freqS_Ts_temp[$freq_AS[$i]]++;}  else {$freqS_Tv_temp[$freq_AS[$i]]++;}}}
						if (('TS' ne $most_common_codon[$i]) and ('TS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if ($nuc_NNN_codon eq 'C')			{$freqS_Ts_temp[$freq_TS[$i]]++;}  else {$freqS_Tv_temp[$freq_TS[$i]]++;}}
																						else									{ if ($most_common_codon[$i] =~ /C/g)	{$freqS_Ts_temp[$freq_TS[$i]]++;}  else {$freqS_Tv_temp[$freq_TS[$i]]++;}}}
						if (('CS' ne $most_common_codon[$i]) and ('CS' ne $posOG[$i])) {if ($most_common_codon[$i] eq 'NNN')	{ if ($nuc_NNN_codon eq 'T')			{$freqS_Ts_temp[$freq_CS[$i]]++;}  else {$freqS_Tv_temp[$freq_CS[$i]]++;}}
																						else									{ if ($most_common_codon[$i] =~ /T/g)	{$freqS_Ts_temp[$freq_CS[$i]]++;}  else {$freqS_Tv_temp[$freq_CS[$i]]++;}}}
						if ('NNN' ne $most_common_codon[$i])	{if ($most_common_codon[$i]=~ /S/g)	{if ($most_common_codon[$i]=~ /G/g) { if ($nuc_NNN_codon eq 'A')	{$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}}
																									 if ($most_common_codon[$i]=~ /A/g) { if ($nuc_NNN_codon eq 'G')	{$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}}
																									 if ($most_common_codon[$i]=~ /C/g) { if ($nuc_NNN_codon eq 'T')	{$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}}
																									 if ($most_common_codon[$i]=~ /T/g) { if ($nuc_NNN_codon eq 'C')	{$freqS_Ts_temp[$freq_NNN[$i]]++;} else {$freqS_Tv_temp[$freq_NNN[$i]]++;}}	}
																}
			
						}	
					 }
				} # elsif outgroup has no NNN
			else {die, "mutation not assigned";}	

	print OUT_DIFF "polytable Repl:\t\t", join ("-", @poly_freq_Rep_temp), "\n";
	print OUT_DIFF "polytable Syn:\t\t", join ("-", @poly_freq_Syn_temp), "\n";
	print OUT_DIFF "polytable U -> U:\t", join ("-", @freq_U_U_temp), "\n";
	print OUT_DIFF "polytable P -> P:\t", join ("-", @freq_P_P_temp), "\n";
	print OUT_DIFF "polytable U -> P:\t", join ("-", @freq_U_P_temp), "\n";
	print OUT_DIFF "polytable P -> U:\t", join ("-", @freq_P_U_temp), "\n";
	print OUT_DIFF "Syn AT -> GC:\t\t", join ("-", @freqS_AT_GC_temp), "\n";
	print OUT_DIFF "Syn AT -> AT:\t\t", join ("-", @freqS_AT_AT_temp), "\n";
	print OUT_DIFF "Syn GC -> GC:\t\t", join ("-", @freqS_GC_GC_temp), "\n";
	print OUT_DIFF "Syn GC -> AT:\t\t", join ("-", @freqS_GC_AT_temp), "\n";
	print OUT_DIFF "Rep AT -> GC:\t\t", join ("-", @freqR_AT_GC_temp), "\n";
	print OUT_DIFF "Rep AT -> AT:\t\t", join ("-", @freqR_AT_AT_temp), "\n";
	print OUT_DIFF "Rep GC -> GC:\t\t", join ("-", @freqR_GC_GC_temp), "\n";
	print OUT_DIFF "Rep GC -> AT:\t\t", join ("-", @freqR_GC_AT_temp), "\n";
	print OUT_DIFF "Syn Ts:\t\t", join ("-", @freqS_Ts_temp), "\n";
	print OUT_DIFF "Syn Tv:\t\t", join ("-", @freqS_Tv_temp), "\n";	

     } # end of going through the three positons..

#	print "polytable  Synonymous min: ", join ("-", @path_min_Syn), "\n";
#	print "polytable Replacement min: ", join ("-", @path_min_Rep), "\n";

# compare path to previous -- if shorter, take it

	$no_polyS_min=$no_polyR_min=0;
	$no_polyS=$no_polyR=0;
	for ($ind=1; $ind<$numseqs; $ind++)
		{
		$no_polyS=$no_polyS + $poly_freq_Syn_temp[$ind];
		$no_polyR=$no_polyR + $poly_freq_Rep_temp[$ind];
		$no_polyS_min=$no_polyS_min+$path_min_Syn[$ind];
		$no_polyR_min=$no_polyR_min+$path_min_Rep[$ind];
		}
	$no_divS=$poly_freq_Syn_temp[$numseqs];
	$no_divR=$poly_freq_Rep_temp[$numseqs];
	$no_divS_min=$path_min_Syn[$numseqs];
	$no_divR_min=$path_min_Rep[$numseqs];
	$no_mut_current= $no_divS + $no_divR + $no_polyS + $no_polyR;
	$no_mut_min= $no_divS_min + $no_divR_min + $no_polyS_min + $no_polyR_min;
	print OUT_DIFF "Synonymous: poly: $no_polyS, divergence: $no_divS\n";
	print OUT_DIFF "Replacement: poly: $no_polyR, divergence: $no_divR\n";
#	print "Synonymous min: poly: $no_polyS_min, divergence: $no_divS_min\n";
#	print "Replacement min: poly: $no_polyR_min, divergence: $no_divR_min\n";

 $identical=1;
 for ($ind=0; $ind<$numseqs; $ind++)
    {
	if (($poly_freq_Syn_temp[$ind] == $path_min_Syn[$ind]) and  ($poly_freq_Rep_temp[$ind] == $path_min_Rep[$ind]))
		{}
	else {$identical=0;}	
	}

  if ($identical==1)
	{}
  elsif ($no_mut_current < $no_mut_min)
	{
	@path_min_Rep = @poly_freq_Rep_temp;
	@path_min_Syn = @poly_freq_Syn_temp;
	@path_min_U_U = @freq_U_U_temp;			@path_min_P_P = @freq_P_P_temp;			@path_min_U_P = @freq_U_P_temp;			@path_min_P_U = @freq_P_U_temp;
	@path_minS_AT_GC = @freqS_AT_GC_temp;	@path_minS_AT_AT = @freqS_AT_AT_temp;	@path_minS_GC_GC = @freqS_GC_GC_temp;	@path_minS_GC_AT = @freqS_GC_AT_temp;
	@path_minR_AT_GC = @freqR_AT_GC_temp;	@path_minR_AT_AT = @freqR_AT_AT_temp;	@path_minR_GC_GC = @freqR_GC_GC_temp;	@path_minR_GC_AT = @freqR_GC_AT_temp;
	@path_minS_Ts = @freqS_Ts_temp;		@path_minS_Tv = @freqS_Tv_temp;
	}
  elsif ($no_divR < $no_divR_min)
	{
	@path_min_Rep = @poly_freq_Rep_temp;
	@path_min_Syn = @poly_freq_Syn_temp;
	@path_min_U_U = @freq_U_U_temp;			@path_min_P_P = @freq_P_P_temp;			@path_min_U_P = @freq_U_P_temp;			@path_min_P_U = @freq_P_U_temp;
	@path_minS_AT_GC = @freqS_AT_GC_temp;	@path_minS_AT_AT = @freqS_AT_AT_temp;	@path_minS_GC_GC = @freqS_GC_GC_temp;	@path_minS_GC_AT = @freqS_GC_AT_temp;
	@path_minR_AT_GC = @freqR_AT_GC_temp;	@path_minR_AT_AT = @freqR_AT_AT_temp;	@path_minR_GC_GC = @freqR_GC_GC_temp;	@path_minR_GC_AT = @freqR_GC_AT_temp;
	@path_minS_Ts = @freqS_Ts_temp;		@path_minS_Tv = @freqS_Tv_temp;
	}
  elsif (($no_divR == $no_divR_min) and ($no_polyR < $no_polyR_min))	
	{
	@path_min_Rep = @poly_freq_Rep_temp;
	@path_min_Syn = @poly_freq_Syn_temp;
	@path_min_U_U = @freq_U_U_temp;			@path_min_P_P = @freq_P_P_temp;			@path_min_U_P = @freq_U_P_temp;			@path_min_P_U = @freq_P_U_temp;
	@path_minS_AT_GC = @freqS_AT_GC_temp;	@path_minS_AT_AT = @freqS_AT_AT_temp;	@path_minS_GC_GC = @freqS_GC_GC_temp;	@path_minS_GC_AT = @freqS_GC_AT_temp;
	@path_minR_AT_GC = @freqR_AT_GC_temp;	@path_minR_AT_AT = @freqR_AT_AT_temp;	@path_minR_GC_GC = @freqR_GC_GC_temp;	@path_minR_GC_AT = @freqR_GC_AT_temp;
	@path_minS_Ts = @freqS_Ts_temp;		@path_minS_Tv = @freqS_Tv_temp;
	}
  elsif ($no_divS < $no_divS_min)
  	{
	@path_min_Rep = @poly_freq_Rep_temp;
	@path_min_Syn = @poly_freq_Syn_temp;
	@path_min_U_U = @freq_U_U_temp;			@path_min_P_P = @freq_P_P_temp;			@path_min_U_P = @freq_U_P_temp;			@path_min_P_U = @freq_P_U_temp;
	@path_minS_AT_GC = @freqS_AT_GC_temp;	@path_minS_AT_AT = @freqS_AT_AT_temp;	@path_minS_GC_GC = @freqS_GC_GC_temp;	@path_minS_GC_AT = @freqS_GC_AT_temp;
	@path_minR_AT_GC = @freqR_AT_GC_temp;	@path_minR_AT_AT = @freqR_AT_AT_temp;	@path_minR_GC_GC = @freqR_GC_GC_temp;	@path_minR_GC_AT = @freqR_GC_AT_temp;
	@path_minS_Ts = @freqS_Ts_temp;		@path_minS_Tv = @freqS_Tv_temp;
	}
  elsif (($no_divS == $no_divS_min) and ($no_polyS < $no_polyS_min))	
	{
	@path_min_Rep = @poly_freq_Rep_temp;
	@path_min_Syn = @poly_freq_Syn_temp;
	@path_min_U_U = @freq_U_U_temp;			@path_min_P_P = @freq_P_P_temp;			@path_min_U_P = @freq_U_P_temp;			@path_min_P_U = @freq_P_U_temp;
	@path_minS_AT_GC = @freqS_AT_GC_temp;	@path_minS_AT_AT = @freqS_AT_AT_temp;	@path_minS_GC_GC = @freqS_GC_GC_temp;	@path_minS_GC_AT = @freqS_GC_AT_temp;
	@path_minR_AT_GC = @freqR_AT_GC_temp;	@path_minR_AT_AT = @freqR_AT_AT_temp;	@path_minR_GC_GC = @freqR_GC_GC_temp;	@path_minR_GC_AT = @freqR_GC_AT_temp;
	@path_minS_Ts = @freqS_Ts_temp;		@path_minS_Tv = @freqS_Tv_temp;
	}

  
 } # end to loop through all possible paths 

	print OUT_DIFF "polytable  Synonymous min: ", join ("-", @path_min_Syn), "\n";
	print OUT_DIFF "polytable Replacement min: ", join ("-", @path_min_Rep), "\n";
	print OUT_DIFF "polytable U -> U min: ", join ("-", @path_min_U_U), "\n";
	print OUT_DIFF "polytable P -> P min: ", join ("-", @path_min_P_P), "\n";
	print OUT_DIFF "polytable U -> P min: ", join ("-", @path_min_U_P), "\n";
	print OUT_DIFF "polytable P -> U min: ", join ("-", @path_min_P_U), "\n";
	print OUT_DIFF "polytable Ts min: ", join ("-", @path_minS_Ts), "\n";
	print OUT_DIFF "polytable Tv min: ", join ("-", @path_minS_Tv), "\n";


# add shortest path to @poly_freq_Syn and @poly_freq_Syn
for ($ind=0; $ind<($numseqs+1); $ind++)
	{
	$poly_freq_Syn[$ind]=$poly_freq_Syn[$ind]+ $path_min_Syn[$ind];
	$poly_freq_Rep[$ind]=$poly_freq_Rep[$ind]+ $path_min_Rep[$ind];
	$freq_U_U[$ind]=$freq_U_U[$ind]+ $path_min_U_U[$ind];
	$freq_P_P[$ind]=$freq_P_P[$ind]+ $path_min_P_P[$ind];
	$freq_U_P[$ind]=$freq_U_P[$ind]+ $path_min_U_P[$ind];
	$freq_P_U[$ind]=$freq_P_U[$ind]+ $path_min_P_U[$ind];	
	$freqS_GC_AT[$ind]=$freqS_GC_AT[$ind]+$path_minS_GC_AT[$ind];
	$freqS_AT_GC[$ind]=$freqS_AT_GC[$ind]+$path_minS_AT_GC[$ind];
	$freqS_GC_GC[$ind]=$freqS_GC_GC[$ind]+$path_minS_GC_GC[$ind];
	$freqS_AT_AT[$ind]=$freqS_AT_AT[$ind]+$path_minS_AT_AT[$ind];
	$freqR_GC_AT[$ind]=$freqR_GC_AT[$ind]+$path_minR_GC_AT[$ind];
	$freqR_AT_GC[$ind]=$freqR_AT_GC[$ind]+$path_minR_AT_GC[$ind];
	$freqR_GC_GC[$ind]=$freqR_GC_GC[$ind]+$path_minR_GC_GC[$ind];
	$freqR_AT_AT[$ind]=$freqR_AT_AT[$ind]+$path_minR_AT_AT[$ind];
	$freqS_Ts[$ind]=$freqS_Ts[$ind]+$path_minS_Ts[$ind];
	$freqS_Tv[$ind]=$freqS_Tv[$ind]+$path_minS_Tv[$ind];
	}
#print OUT_DIFF "polytable Synonymous: ",	join ("-", @poly_freq_Syn), "\n";
#print OUT_DIFF "polytable Replacement: ",	join ("-", @poly_freq_Rep), "\n";
	  	  
} # if $unique_codons>2
        
	   
}  # loop for no_GAP --- if ($switch_gap==0)	
}  # loop for all codons $pos


$no_polyS=$no_polyR=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$no_polyS=$no_polyS+$poly_freq_Syn[$ind];
	$no_polyR=$no_polyR+$poly_freq_Rep[$ind];
	}
$no_divS=$poly_freq_Syn[$numseqs];
$no_divR=$poly_freq_Rep[$numseqs];

$no_divP_U=$freq_P_U[$numseqs];		$no_divU_P=$freq_U_P[$numseqs];		$no_divP_P=$freq_P_P[$numseqs];		$no_divU_U=$freq_U_U[$numseqs];

$no_divS_AT_GC=$freqS_AT_GC[$numseqs];		$no_divS_GC_AT=$freqS_GC_AT[$numseqs];		$no_divS_AT_AT=$freqS_AT_AT[$numseqs];		$no_divS_GC_GC=$freqS_GC_GC[$numseqs];
$no_divR_AT_GC=$freqR_AT_GC[$numseqs];		$no_divR_GC_AT=$freqR_GC_AT[$numseqs];		$no_divR_AT_AT=$freqR_AT_AT[$numseqs];		$no_divR_GC_GC=$freqR_GC_GC[$numseqs];

$no_divS_Ts=$freqS_Ts[$numseqs];	$no_divS_Tv=$freqS_Tv[$numseqs];

if ($freq_cut_off != 0)
	{
	$no_polyS_freq=$no_polyR_freq=0;
	$no_polyP_U_freq= $no_polyU_P_freq= $no_polyP_P_freq= $no_polyU_U_freq=0;
	
	for ($ind=1; $ind<$numseqs-1; $ind++)
		{
		if ( ($ind/($numseqs-1)) > $freq_cut_off)
			{
			$no_polyS_freq=$no_polyS_freq+$poly_freq_Syn[$ind];
			$no_polyR_freq=$no_polyR_freq+$poly_freq_Rep[$ind];
			
			$no_polyP_U_freq=$no_polyP_U_freq+$freq_P_U[$ind];
			$no_polyU_P_freq=$no_polyU_P_freq+$freq_U_P[$ind];
			$no_polyP_P_freq=$no_polyP_P_freq+$freq_P_P[$ind];	
			$no_polyU_U_freq=$no_polyU_U_freq+$freq_U_U[$ind];
			}
		}
	}

$no_syn_codons=$no_syn_codons/$numseqs;
$no_rep_codons=$no_rep_codons/$numseqs;
$no_syn_fourfold_codons=$no_syn_fourfold_codons/$numseqs;

# calculate pi
$pi_syn_total=$pi_rep_total=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$pi_syn[$ind]=(2*($ind/($numseqs-1))*(1-($ind/($numseqs-1))))*$poly_freq_Syn[$ind];
	$pi_rep[$ind]=(2*($ind/($numseqs-1))*(1-($ind/($numseqs-1))))*$poly_freq_Rep[$ind];
	$pi_syn_total=$pi_syn_total+$pi_syn[$ind];
	$pi_rep_total=$pi_rep_total+$pi_rep[$ind];
	}
$pi_syn_total=$pi_syn_total*(($numseqs-1)/($numseqs-2));	
$pi_rep_total=$pi_rep_total*(($numseqs-1)/($numseqs-2));
#print "no_syn_codons ", $no_syn_codons;
$pi_syn_site=$pi_syn_total/$no_syn_codons;
$pi_rep_site=$pi_rep_total/$no_rep_codons;	

$pi_JC_syn= -0.75*log(1-(4/3)*$pi_syn_site);
$pi_JC_rep= -0.75*log(1-(4/3)*$pi_rep_site);


# calculate Dxy
$freq_syn=$freq_rep=0;
$freq_syn_Ts=$freq_syn_Tv=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$freq_syn_temp[$ind]=$ind/($numseqs-1)*$poly_freq_Syn[$ind];
	$freq_rep_temp[$ind]=$ind/($numseqs-1)*$poly_freq_Rep[$ind];
	$freq_syn=$freq_syn+$freq_syn_temp[$ind];
	$freq_rep=$freq_rep+$freq_rep_temp[$ind];
	
	$freq_syn_Ts_temp[$ind]=$ind/($numseqs-1)*$freqS_Ts[$ind];
	$freq_syn_Ts=$freq_syn_Ts+$freq_syn_Ts_temp[$ind];
	$freq_syn_Tv_temp[$ind]=$ind/($numseqs-1)*$freqS_Tv[$ind];
	$freq_syn_Tv=$freq_syn_Tv+$freq_syn_Tv_temp[$ind];
	}
	
$Dxy_syn = ($no_divS+$freq_syn)/$no_syn_codons;
$Dxy_rep = ($no_divR+$freq_rep)/$no_rep_codons;

$Dxy_syn_Ts = ($no_divS_Ts+$freq_syn_Ts)/$no_syn_codons;	# P
$Dxy_syn_Tv = ($no_divS_Tv+$freq_syn_Tv)/$no_syn_codons;	# Q


# $D_fourfold = $no_syn_fourfold_div/$no_syn_fourfold_codons;

$Dxy_JC_syn=  -0.75*log(1-(4/3)*$Dxy_syn);
$Dxy_JC_rep=  -0.75*log(1-(4/3)*$Dxy_rep);
$D_JC_fourfold=  -0.75*log(1-(4/3)*$D_fourfold);

$Dxy_Kimura_syn=  0.5*log(1/(1 - (2*$Dxy_syn_Ts) - $Dxy_syn_Tv)) + (0.25*log(1/(1 - (2*$Dxy_syn_Tv))));

# calculate Fay&Wu's H
$pi_H_syn_total=$pi_H_rep_total=0;
for ($ind=1; $ind<$numseqs-1; $ind++)
	{
	$pi_H_syn[$ind]=2*$ind*$ind*$poly_freq_Syn[$ind];
	$pi_H_syn_total=$pi_H_syn_total+$pi_H_syn[$ind];

	$pi_H_rep[$ind]=2*$ind*$ind*$poly_freq_Rep[$ind];
	$pi_H_rep_total=$pi_H_rep_total+$pi_H_rep[$ind];
	}
$pi_H_syn_total=$pi_H_syn_total*(1/(($numseqs-1)*($numseqs-2)));
$FayWu_H_syn=$pi_syn_total-$pi_H_syn_total;
$pi_H_rep_total=$pi_H_rep_total*(1/(($numseqs-1)*($numseqs-2)));
$FayWu_H_rep=$pi_rep_total-$pi_H_rep_total;

#**************** calculate TajD

$numseqs_ingroup=$numseqs-1;
open INPUT, "./tajd ". $numseqs_ingroup . " " . $no_polyS . " " . $pi_syn_total . " " .$poly_freq_Syn[1] . " |";
# print "./tajd  $numseqs_ingroup  $no_polyS  $pi_syn_total  $poly_freq_Syn[1]\n";
	my $line1 = <INPUT>; chomp $line1;
	my $line2 = <INPUT>; chomp $line2;
	($junk, $TajD_syn)= split(/=/, $line1);
	($junk, $FuLiD_syn)= split(/=/, $line2);
	

open INPUT, "./tajd ". $numseqs_ingroup . " " . $no_polyR . " " . $pi_rep_total . " " .$poly_freq_Rep[1] . " |";
# print "./tajd  $numseqs_ingroup  $no_polyR  $pi_rep_total  $poly_freq_Rep[1]\n";
	my $line1 = <INPUT>; chomp $line1;
	my $line2 = <INPUT>; chomp $line2;
	($junk, $TajD_rep)= split(/=/, $line1);
	($junk, $FuLiD_rep)= split(/=/, $line2);

pop @poly_freq_Syn; pop @poly_freq_Syn;
pop @poly_freq_Rep; pop @poly_freq_Rep;

pop @freq_P_U; pop @freq_P_U;	
pop @freq_U_P; pop @freq_U_P;	
pop @freq_P_P; pop @freq_P_P;	
pop @freq_U_U; pop @freq_U_U;

pop @freqS_AT_GC; pop @freqS_AT_GC; pop @freqR_AT_GC; pop @freqR_AT_GC;
pop @freqS_GC_AT; pop @freqS_GC_AT; pop @freqR_GC_AT; pop @freqR_GC_AT;
pop @freqS_AT_AT; pop @freqS_AT_AT; pop @freqR_AT_AT; pop @freqR_AT_AT;
pop @freqS_GC_GC; pop @freqS_GC_GC; pop @freqR_GC_GC; pop @freqR_GC_GC;

pop @freqS_Ts; pop @freqS_Ts; pop @freqS_Tv; pop @freqS_Tv; 

$poly_freq_Syn[0]=$no_divS;
$poly_freq_Rep[0]=$no_divR;

$freq_P_U[0]=$no_divP_U;	$freq_U_P[0]=$no_divU_P;	$freq_P_P[0]=$no_divP_P;	$freq_U_U[0]=$no_divU_U;
$freqS_AT_GC[0]=$no_divS_AT_GC;		$freqS_GC_AT[0]=$no_divS_GC_AT;		$freqS_AT_AT[0]=$no_divS_AT_AT;		$freqS_GC_GC[0]=$no_divS_GC_GC; 
$freqR_AT_GC[0]=$no_divR_AT_GC;		$freqR_GC_AT[0]=$no_divR_GC_AT;		$freqR_AT_AT[0]=$no_divR_AT_AT;		$freqR_GC_GC[0]=$no_divR_GC_GC; 

$freqS_Ts[0]=$no_divS_Ts;	$freqS_Tv[0]=$no_divS_Tv;

print OUT2 $file, "\t", $sequence_names[0], "\t", $numseqs-1, "\t", $no_syn_codons,"\t", $no_polyS, "\t", $no_polyS_freq, "\t", $no_divS, "\t", $pi_syn_site, "\t", $pi_JC_syn, "\t", $Dxy_syn, "\t", $Dxy_JC_syn, "\t", $Dxy_Kimura_syn, "\t",  $TajD_syn, "\t",  $FuLiD_syn,  "\t", $FayWu_H_syn,  "\t"; 
print OUT2													  "\t", $no_rep_codons,"\t", $no_polyR, "\t", $no_polyR_freq, "\t", $no_divR, "\t", $pi_rep_site, "\t", $pi_JC_rep, "\t", $Dxy_rep, "\t", $Dxy_JC_rep, "\t",						 $TajD_rep, "\t",  $FuLiD_rep,  "\t", $FayWu_H_rep,  "\n";

$GC_three=$GC_three/$third_pos_count;
$GC_three_invariant=$GC_three_invariant/$third_pos_count_invariant;

$FOP = $FOP/($FOP+$FNOP);

# print  "\nlocus ", $file, " numseqs: ", $numseqs, "\n";
# print "Num of Syn Sites: ", $no_syn_codons, "\n";
# print "Num of Rep Sites: ", $no_rep_codons, "\n";
# print  "Synonymous Poly: $no_polyS, Synonymous Divergence: $no_divS \n"; 
# print  "Replacement Poly: $no_polyR, Replacement Divergence: $no_divR \n"; 		
# print "Pairwise Syn diversity: $pi_syn_site \t JC corr:  $pi_JC_syn  \n";
# print "Pairwise Rep diversity: $pi_rep_site \t JC corr:  $pi_JC_rep \n";
# print "Dxy_syn ",  $Dxy_syn, "\tDxy_JC_syn  ", $Dxy_JC_syn , "\n";
# print "Dxy_rep ",  $Dxy_rep, "\tDxy_JC_rep  ", $Dxy_JC_rep , "\n";
# print "Syn Taj D is : $TajD_syn, Fu&Li D is $FuLiD_syn\n";
# print "Rep Taj D is : $TajD_rep, Fu&Li D is $FuLiD_rep\n";
# print  "Synonymous Poly above a frequency of $freq_cut_off: $no_polyS_freq  \n"; 
# print  "Replacement Polyabove a frequency of $freq_cut_off: $no_polyR_freq  \n"; 		

print OUT $file, "_Syn\t", join ("\t", @poly_freq_Syn), "\n";
print OUT $file, "_Rep\t" , join ("\t", @poly_freq_Rep), "\n";

print OUT3 $file, "\t P->U\t", $no_polyP_U_freq, "\t", join ("\t", @freq_P_U), "\n";
print OUT3 $file, "\t U->P\t", $no_polyU_P_freq, "\t", join ("\t", @freq_U_P), "\n";
print OUT3 $file, "\t P->P\t", $no_polyP_P_freq, "\t", join ("\t", @freq_P_P), "\n";
print OUT3 $file, "\t U->U\t", $no_polyU_U_freq, "\t", join ("\t", @freq_U_U), "\n";


print OUT4 $file, "\t %GC3(all_sites) \t", $GC_three, "\n";
print OUT4 $file, "\t %GC3(invariant) \t", $GC_three_invariant, "\n";
print OUT4 $file, "\t FOP \t", $FOP, "\n";
print OUT4 $file, "\t Syn_AT->GC \t", join ("\t", @freqS_AT_GC), "\n";
print OUT4 $file, "\t Syn_GC->AT \t", join ("\t", @freqS_GC_AT), "\n";
print OUT4 $file, "\t Syn_AT->AT \t", join ("\t", @freqS_AT_AT), "\n";
print OUT4 $file, "\t Syn_GC->GC \t", join ("\t", @freqS_GC_GC), "\n";
print OUT4 $file, "\t Rep_AT->GC \t", join ("\t", @freqR_AT_GC), "\n";
print OUT4 $file, "\t Rep_GC->AT \t", join ("\t", @freqR_GC_AT), "\n";
print OUT4 $file, "\t Rep_AT->AT \t", join ("\t", @freqR_AT_AT), "\n";
print OUT4 $file, "\t Rep_GC->GC \t", join ("\t", @freqR_GC_GC), "\n";

print OUT5 $file, "\t Syn_Ts \t", join ("\t", @freqS_Ts), "\n";
print OUT5 $file, "\t Syn_Tv \t", join ("\t", @freqS_Tv), "\n";

} # loop foreach file







# *************************************************************************
# SUBROUTINE codon_processor
# 
# input is an array (@codon_list_pop) that returns a modified codon_list


 sub codon_processor {

 (@codon_list_pop) = @_;
 $complexity = pop @codon_list_pop;

# print what came in
# print join("\t", @codon_list_pop), "\n";
# print "complex codon: $complexity \n";

# make a list @codon_list that contains only unique codons
# use each element in this list as the outgroup codons - if the outgroup codon is not present in the population -
# and go through the entire program 
@codon_list =();   # contains only unique codons
@temp_codon_list = sort @codon_list_pop;
push @codon_list, $temp_codon_list[0];
for ($x=0; $x<(scalar(@temp_codon_list)-1); $x++)
	   {
		if ($temp_codon_list[$x] ne $temp_codon_list[$x+1])
			{
			push @codon_list, $temp_codon_list[$x+1];
			}
		}
# print join("\t", @codon_list), "\n";

#*** check whether outgroup codon is found in ingroup***
# if yes, use outgroup sequence as ancestral codon and only loop through SUBROUTINE codon_processor once
# if not, go through subroutine multiple times, and take path that minimises the number of replacement
# mutations in ingroup
 

# if ougroup codon either is found in ingroup, or it is a fixed difference and no polymorphism -- go through loop only once
if ($complexity eq 'nein')
   {
	$switch_subroutine=1;
	}	
else 
	{            
	$switch_subroutine = scalar(@codon_list);  # divergent site but polymorphic in ingroup
	}

# print "no of unique codons: ", scalar(@codon_list), "\n";
# print "switch = $switch_subroutine, $complexity \n";
# print "\n";

#my @codon = ();
#my @codon_list = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]); 
$out=0;
@final_hashes = ();



for ($out=0; $out<$switch_subroutine; $out++)
  {
# if we go through loop multiple times, $codon_OG is the $out codon in the unique codon list
# if we go through only once, the outgroup $codon_list_pop[0] is the $codon_OG 
if ($switch_subroutine == 1)
   {
	 $codon_OG = $codon_list_pop[0];
	}	
else 
	{            
	$codon_OG = $codon_list[$out];  # codon not present
	}
#  print "outgroup used: ", $codon_OG, "\n";



$x=$y=0; 
$i=$j=$k=0;
$bases=4;
$no_codons=64;

# define mut = all 1 2 and 3 steps away from outgrup codon "000"
   
	@mut = ('000', '300', '200', '100',  '030', '020', '010', '003', '002', '001', 	 
			'300', '330', '320', '310', '303', '302', '301',
			'200', '230', '220', '210', '203', '202', '201',
			'100', '130', '120', '110', '103', '102', '101',
			'030', '330', '230', '130', '033', '032', '031', 
			'020', '320', '220', '120', '023', '022', '021',
			'010', '310', '210', '110', '013', '012', '011',
			'003', '303', '203', '103', '033', '023', '013',
			'002', '302', '202', '102', '032', '022', '012',
			'001', '301', '201', '101', '031', '021', '011',
			'330', '333', '332', '331',
			'320', '323', '322', '321',
			'310', '313', '312', '311',
			'303', '333', '323', '313',
			'302', '332', '322', '312',
			'301', '331', '321', '311',
			'230', '233', '232', '231',
			'220', '223', '222', '221',
			'210', '213', '212', '211',
			'203', '233', '223', '213',
			'202', '232', '222', '212',
			'201', '231', '221', '211',	
			'130', '133', '132', '131',
			'120', '123', '122', '121',
			'110', '113', '112', '111',
			'103', '133', '123', '113',
			'102', '132', '122', '112',
			'101', '131', '121', '111',
			'033', '333', '233', '133',
			'032', '332', '232', '132',
			'031', '331', '231', '131',
			'023', '323', '223', '123',
			'022', '322', '222', '122',
			'021', '321', '221', '121',
			'013', '313', '213', '113',
			'012', '312', '212', '112',
			'011', '311', '211', '111');
			   
			   
$array_length = scalar(@mut);
# print "array length: ", $array_length , "\n";

#########################################################################################
# assign bases assuming that outgroup has base "000"

for ($x=0; $x<3; $x++)
   {
    
	for ($y=0; $y<$array_length; $y++)
		{
			if (substr($codon_OG,$x,1) eq 'A')				# assign A ancestral
			{
			if (substr($mut[$y], $x, 1) eq '0')
				{substr($mut[$y], $x, 1) = 'A';}
			if (substr($mut[$y], $x, 1) eq '1')
				{substr($mut[$y], $x, 1) = 'T';}
			if (substr($mut[$y], $x, 1) eq '2')
				{substr($mut[$y], $x, 1) = 'G';}
			if (substr($mut[$y], $x, 1) eq '3')
				{substr($mut[$y], $x, 1) = 'C';}
			}
			
		if (substr($codon_OG,$x,1) eq 'G')				# assign G ancestral
			{
	   	 
			if (substr($mut[$y], $x, 1) eq '0')
				{substr($mut[$y], $x, 1) = 'G';}
			if (substr($mut[$y], $x, 1) eq '1')
				{substr($mut[$y], $x, 1) = 'T'}
			if (substr($mut[$y], $x, 1) eq '2')
				{substr($mut[$y], $x, 1) = 'A';}
			if (substr($mut[$y], $x, 1) eq '3')
				{substr($mut[$y], $x, 1) = 'C';}
			}
			
			if (substr($codon_OG,$x,1) eq 'T')				# assign T ancestral
			{
	   	 
			if (substr($mut[$y], $x, 1) eq '0')
				{substr($mut[$y], $x, 1) = 'T';}
			if (substr($mut[$y], $x, 1) eq '1')
				{substr($mut[$y], $x, 1) = 'G';}
			if (substr($mut[$y], $x, 1) eq '2')
				{substr($mut[$y], $x, 1) = 'A';}
			if (substr($mut[$y], $x, 1) eq '3')
				{substr($mut[$y], $x, 1) = 'C';}
			}
			
			if (substr($codon_OG,$x,1) eq 'C')				# assign C ancestral
			{
	   	 
			if (substr($mut[$y], $x, 1) eq '0')
				{substr($mut[$y], $x, 1) = 'C';}
			if (substr($mut[$y], $x, 1) eq '1')
				{substr($mut[$y], $x, 1) = 'G';}
			if (substr($mut[$y], $x, 1) eq '2')
				{substr($mut[$y], $x, 1) = 'A';}
			if (substr($mut[$y], $x, 1) eq '3')
				{substr($mut[$y], $x, 1) = 'T';}
			}

		 }
	}

#	print @mut, "\n\n";
######################################################################################################################
#
# make hashes that are one mutational step away (%muthash1), two mutational steps away (%muthash2) or three mutational
# steps away (%muthash3) from the outgroup sequence "NNN"
#


%muthash1 = ($mut[0] => [$mut[1], $mut[2], $mut[3], $mut[4], $mut[5], $mut[6], $mut[7], $mut[8], $mut[9]]);

%muthash2 = ($mut[10] => [$mut[11], $mut[12], $mut[13], $mut[14], $mut[15], $mut[16]],
             $mut[17] => [$mut[18], $mut[19], $mut[20], $mut[21], $mut[22], $mut[23]],
			 $mut[24] => [$mut[25], $mut[26], $mut[27], $mut[28], $mut[29], $mut[30]],
			 $mut[31] => [$mut[32], $mut[33], $mut[34], $mut[35], $mut[36], $mut[37]],
			 $mut[38] => [$mut[39], $mut[40], $mut[41], $mut[42], $mut[43], $mut[44]],
			 $mut[45] => [$mut[46], $mut[47], $mut[48], $mut[49], $mut[50], $mut[51]],
			 $mut[52] => [$mut[53], $mut[54], $mut[55], $mut[56], $mut[57], $mut[58]], 
             $mut[59] => [$mut[60], $mut[61], $mut[62], $mut[63], $mut[64], $mut[65]], 
             $mut[66] => [$mut[67], $mut[68], $mut[69], $mut[70], $mut[71], $mut[72]] );
			 
%muthash3 = ($mut[73] => [$mut[74], $mut[75], $mut[76]],
             $mut[77] => [$mut[78], $mut[79], $mut[80]],
			 $mut[81] => [$mut[82], $mut[83], $mut[84]],
			 $mut[85] => [$mut[86], $mut[87], $mut[88]],
			 $mut[89] => [$mut[90], $mut[91], $mut[92]],
			 $mut[93] => [$mut[94], $mut[95], $mut[96]],
			 $mut[97] => [$mut[98], $mut[99], $mut[100]], 
             $mut[101] => [$mut[102], $mut[103], $mut[104]], 
             $mut[105] => [$mut[106], $mut[107], $mut[108]], 
             $mut[109] => [$mut[110], $mut[111], $mut[112]], 
             $mut[113] => [$mut[114], $mut[115], $mut[116]], 
             $mut[117] => [$mut[118], $mut[119], $mut[120]], 
             $mut[121] => [$mut[122], $mut[123], $mut[124]], 
             $mut[125] => [$mut[126], $mut[127], $mut[128]], 
             $mut[129] => [$mut[130], $mut[131], $mut[132]], 
             $mut[133] => [$mut[134], $mut[135], $mut[136]], 
             $mut[137] => [$mut[138], $mut[139], $mut[140]], 
             $mut[141] => [$mut[142], $mut[143], $mut[144]], 
             $mut[145] => [$mut[146], $mut[147], $mut[148]], 
             $mut[149] => [$mut[150], $mut[151], $mut[152]], 
             $mut[153] => [$mut[154], $mut[155], $mut[156]], 
             $mut[157] => [$mut[158], $mut[159], $mut[160]], 
             $mut[161] => [$mut[162], $mut[163], $mut[164]], 
             $mut[165] => [$mut[166], $mut[167], $mut[168]], 
             $mut[169] => [$mut[170], $mut[171], $mut[172]], 
             $mut[173] => [$mut[174], $mut[175], $mut[176]], 
             $mut[177] => [$mut[178], $mut[179], $mut[180]] 
			 );

# print "key: ", keys %muthash1, "\n";
# foreach $item (sort (keys(%muthash1)))
# { print "Key $item: ", join ("--", @{$muthash1{$item}}), "\n";}
# foreach $item (sort (keys(%muthash2)))
#  {print "Key $item: ", join ("--", @{$muthash2{$item}}), "\n";}
# foreach $item (sort (keys(%muthash3)))
# {print "Key $item: ", join ("--", @{$muthash3{$item}}), "\n";}


$length_muthash1 = 9;
$length_muthash2 = 54;
$length_muthash3 = 81;

##############################################################################################
# go through all codons in the population and assign them to different mutational classes 
# (i.e. 1, 2 or 3 mutational steps away), and put the codons that are different from the outgroup
# in temporary hashes
##################

$length_codonlist = scalar(@codon_list);  #codon list is entered via the commandline for now
# print "length_codonlist: ", $length_codonlist, "\n";

my @paths3 = ();  # these arrays store possible pathways from a codon to the ancestral codon
my @paths2 = (); 
my @paths1 = ();  

my %mut1_temp = ();  # these arrays store extant codons in the codon_list
my %mut2_temp = ();
my %mut3_temp = ();


$num_paths1=$num_paths2=$num_paths3 = 0;  # keeps track to total number of paths

$y=0;	# loop through all codons in codon_list (from codon 1 to codon n  ----   changed this to zero!!!!)
for ($y=0; $y< $length_codonlist; $y++)
{
	$x=0;
	$switch='no';
	for ($x=0; $x<$length_muthash1; $x++)   # check for that codon in muthash1
		{
		# print "test: ", $muthash1{$codon_OG}[$x], "\n";
	
	if ($muthash1{$codon_OG}[$x] eq $codon_list[$y])
			{
			$switch='yes';
			# assign to a hash, the value of which is always the ancestral state
			$mut1_temp{$codon_list[$y]} = $codon_OG; 			
			}
		 }		 
		# print "is ", $codon_list[$y], " an element of muthash1?: ", $switch, "\n"; 

	if ($switch eq 'no')  # if that codon was not in muthash1 go through muthash2
			{    
			$x=0;
			$i=0;
			for ($x=0; $x<$length_muthash2; $x++)   # search through muthash2
				{
              # print "test: ", $muthash2{$codon_OG}[$x], "\n";
			
				foreach $item (sort(keys(%muthash2)))
					{
					if ($muthash2{$item}[$x] eq $codon_list[$y])
						{
						$switch='yes';
						# since you found it, now assign it to a hash, the value of which is the ancestral state
						# find key1 and key2
						$key_temp[$i] = $item;
				
						#  print "item: ", $i, "\t", $item, "\n";
				   
						$i++;
						$mut2_temp{$codon_list[$y]} = [$key_temp[0],$key_temp[1]];
						}
					}
                }
			# print "is ", $codon_list[$y], " an element of muthash2?: ", $switch, "\n"; 
			}
		
	if ($switch eq 'no')   # if that codon was not in muthash1 or 2 go through muthash3
			{    
			$x=0;
			$i=0;
			for ($x=0; $x<$length_muthash3; $x++)   # search through muthash3
              {
              # print "test: ", $muthash3{$codon_OG}[$x], "\n";
	
			foreach $item (sort(keys(%muthash3)))
				{
				if ($muthash3{$item}[$x] eq $codon_list[$y])
	               {
				   $switch='yes';
				   # find key1, key2 and key3
				   $key_temp[$i] = $item;
				   
				   # print "item: ", $i, "\t", $item, "\n";
				   $i++;
				   $mut3_temp{$codon_list[$y]} = [$key_temp[0],$key_temp[1], $key_temp[2]];
	               }
				}
               }
			   # print "is ", $codon_list[$y], " an element of muthash3?: ", $switch, "\n"; 
			}
	
	
	
}

#	print "\nall extant codons and their ancestors:\n";
#   print the keys in mut1_temp
		@allkeys_tempmut1 = keys(%mut1_temp);
		foreach $key (@allkeys_tempmut1)
		{
#		print "key is ", $key, " and members of the key are:\t", $mut1_temp{$key}, "\n";
		}

	# print all the keys in mut2_temp
		@allkeys_tempmut2 = keys(%mut2_temp);
		foreach $key (@allkeys_tempmut2)
		{
#		print "key is ", $key, " and members of the key are:\t", $mut2_temp{$key}[0], "\t", $mut2_temp{$key}[1], "\n";
		}	
				
	# print all the keys in mut3_temp
		@allkeys_tempmut3 = keys(%mut3_temp);
		foreach $key (@allkeys_tempmut3)
		{
#		print "key is ", $key, " and members of the key are:\t", $mut3_temp{$key}[0], "\t", $mut3_temp{$key}[1], "\t", $mut3_temp{$key}[2],"\n";
		}

#************************************************************
# assigning all of the codons to mut_temp arrays is complete
#*************************************************************
#**************************************************************************************************************
# all extant codons are saved in %mut1_temp, %mut2_temp and %mut3_temp - codons are stored as hash keys, 
# and the immediate ancestral codons are saved as array elements
# (i.e. 3 elements for %mut3_temp, 2 elements for %mut2_temp, 1 element for %mut1_temp)
# start with most the complicated codons (i.e. %mut3_temp) and make array to reconstruct ancestry of these codons
# add all extant codons from %mut2_temp and consider only those path that are supported by 

 # -------------------------------------------------------------------
 # task 2 -- find shortest path back to the outgroup
 # make an array "@paths3" that represents possible paths back to the ancestral sequence for the mut3 category codons

foreach $key3 (sort(keys(%mut3_temp)))
{

	$z=0;
	for ($z=0; $z<3; $z++)
	{
	$element = $mut3_temp{$key3}[$z];

	# print "element ", $z, " = ", $element, "\n";

	$x=0;
	$i=0;	
	for ($x=0; $x<$length_muthash2; $x++)   # search through muthash2
		{
		foreach $item (sort(keys(%muthash2)))
				{
				if ($muthash2{$item}[$x] eq $element)
	               {
					# find key1 and key2
				   $key_temp[$z][$i] = $item;
				   $i++;
	               }
				}
			}
	#  print "descendant ", $key3, " and ancestors are:", $key_temp[$z][0], "\t", $key_temp[$z][1], "\n";
	
	}
	push(@paths3 , [$key3, $mut3_temp{$key3}[0], $key_temp[0][0]]);
	push(@paths3 , [$key3, $mut3_temp{$key3}[0], $key_temp[0][1]]);
	push(@paths3 , [$key3, $mut3_temp{$key3}[1], $key_temp[1][0]]);
	push(@paths3 , [$key3, $mut3_temp{$key3}[1], $key_temp[1][1]]);	
	push(@paths3 , [$key3, $mut3_temp{$key3}[2], $key_temp[2][0]]);
	push(@paths3 , [$key3, $mut3_temp{$key3}[2], $key_temp[2][1]]);
	
	$num_paths3+=6;
}

# now make an array "@paths2" that represents possible paths back to the ancestral sequence for the mut2 category codons
foreach $key2 (sort(keys(%mut2_temp)))
	{
	push(@paths2 , [$key2, $mut2_temp{$key2}[0]]);
	push(@paths2 , [$key2, $mut2_temp{$key2}[1]]);
	$num_paths2+=2;
	}	

# now make an array "@paths1" that represents possible paths back to the ancestral sequence for the mut1 category codons
foreach $key1 (sort(keys(%mut1_temp)))
	{
	push(@paths1 , $key1);
	$num_paths1++;
	}	


# ********************************************************************************
# you're done baby - print all the reconstracted paths, if you feel so inclined  *

#	print "\nnumber of 3 step paths ", $num_paths3, "\n";		
#	for ($w=0; $w<$num_paths3; $w++)
#	{
#	$path_temp = $paths3[$w][0]." ".$paths3[$w][1]." ".$paths3[$w][2];
#	print "path ", $w, "\t",  $path_temp, "\n";
#	}
		
#	print "\nnumber of 2 step paths ", $num_paths2, "\n";		
#	for ($w=0; $w<$num_paths2; $w++)
#	{
#	$path_temp = $paths2[$w][0]." ".$paths2[$w][1], "\n";;
#	print "path ", $w, "\t",  $path_temp, "\n";
#	}
	
#	print "\nnumber of 1 step paths ", $num_paths1, "\n";
#	for ($w=0; $w<$num_paths1; $w++)
#	{
#	print "path ", $w, "\t",  $paths1[$w], "\n";
#	}


# ****************************************************************
# task 3 use subroutine combiner make one table for each combination of paths
# count the number of mutational steps
# assign R vs S to those steps

$num_3 = scalar(keys(%mut3_temp));
$num_2 = scalar(keys(%mut2_temp));
$num_1 = scalar(keys(%mut1_temp));

# print "\nnum step 3 codons: = $num_3 , num step 2 codons: = $num_2\n\n";

@tables = combiner($num_3, $num_2);
$num_comp2 = (6 ** $num_3) * (2 ** $num_2);
#for $x (0..($num_comp2-1))
# {
#  for $y (0..($num_3+$num_2-1))
#	{
#	 print $tables[$x][$y] , "\t";
#	}
#  print "\n";
# }


# *************************************************************************
# task 4 calculate the number of mutation steps in each table
# & assign synonymous and nonsynonymous changes
# print path with the fewest number of changes and maximum number
# of synonymous changes to @shortest_path 
# *************************************************************************

$columns_table = scalar @{@tables[0]};
@shortest_path = ();
$min_length = $num_3*3+$num_2*2+$num_1;
# print "min_length: $min_length\n";
@reconstr_codons = ();
$maxS=$countS=0;

# THE BIG LOOP ***************************************************************
for $i (0..($num_comp2-1))  # loop through all possible combinations
{

%newlist = ();
$newlist{$codon_OG} = (0);
@mut = ();
$z=0;

# go through all 3-step codons and reconstruct ancestor
for $which_codon3 (0..($num_3-1)) # 3-step mutations
   {
   $reconstr ='NNNNNNNNN';
   $which_path3 = $tables[$i][$which_codon3];
   $which_path3 = $which_path3 + $which_codon3*6;

# reconstruct a 3-step codon ($paths3[$which_path3][0]) from the 2-step codon in its path ($paths3[$which_path3][1]) 
	   $aa_anc[$i]=codon2aa($paths3[$which_path3][1]);
	   $codon_anc[$i]=codonbias($paths3[$which_path3][1]);
	   $aa_des[$i]=codon2aa($paths3[$which_path3][0]);
	   $codon_des[$i]=codonbias($paths3[$which_path3][0]);
	   if ($aa_anc[$i] eq $aa_des[$i])
	         {
			 if   (($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'P')) {$state ='SP';}
			 elsif(($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'U')) {$state ='SD';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'P')) {$state ='SB';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'U')) {$state ='SU';} 
			 else  {$state ='S_';}
			}
	   else  {$state = 'R_';}

	   if (substr(($paths3[$which_path3][0]), 0, 1) ne substr(($paths3[$which_path3][1]), 0, 1))
		 {
		 $reconstr = substr(($paths3[$which_path3][0]), 0, 1).$state.substr($reconstr, 3, 6);
		 }
	   if (substr(($paths3[$which_path3][0]), 1, 1) ne substr(($paths3[$which_path3][1]), 1, 1))
		 {
		 $reconstr = substr($reconstr, 0, 3).substr(($paths3[$which_path3][0]), 1, 1).$state.substr($reconstr, 6, 3);
		 }
	   if (substr(($paths3[$which_path3][0]), 2, 1) ne substr(($paths3[$which_path3][1]), 2, 1))
		 {
		 $reconstr = substr($reconstr, 0, 6).substr(($paths3[$which_path3][0]), 2, 1).$state;
		 }

   if (exists $newlist{$paths3[$which_path3][0]}) {}
   else{
	   $newlist{$paths3[$which_path3][0]} = 0;
	   push @mut, $state;	  
	   #print "$paths3[$which_path3][0] (aa $aa_anc[$i]) to $paths3[$which_path3][1] (aa $aa_des[$i]) is a ", $state, " mut!\n";	  
	   }
	   
# reconstruct the {reconstructed} 2-step codon ($paths3[$which_path3][1]) from the 1-step codon in its path ($paths3[$which_path3][2]) 
	   $aa_anc[$i]=codon2aa($paths3[$which_path3][2]);
	   $codon_anc[$i]=codonbias($paths3[$which_path3][2]);
	   $aa_des[$i]=codon2aa($paths3[$which_path3][1]);
	   $codon_des[$i]=codonbias($paths3[$which_path3][1]);
	   if ($aa_anc[$i] eq $aa_des[$i])
	         {
			 if   (($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'P')) {$state ='SP';}
			 elsif(($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'U')) {$state ='SD';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'P')) {$state ='SB';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'U')) {$state ='SU';} 
			 else  {$state ='S_';}
			}
	   else  {$state = 'R_';}
	   if (substr(($paths3[$which_path3][1]), 0, 1) ne substr(($paths3[$which_path3][2]), 0, 1))
		 {
		 $reconstr = substr(($paths3[$which_path3][1]), 0, 1).$state.substr($reconstr, 3, 6);
		 }
	   if (substr(($paths3[$which_path3][1]), 1, 1) ne substr(($paths3[$which_path3][2]), 1, 1))
		 {
		 $reconstr = substr($reconstr, 0, 3).substr(($paths3[$which_path3][1]), 1, 1).$state.substr($reconstr, 6, 3);
		 }
	   if (substr(($paths3[$which_path3][1]), 2, 1) ne substr(($paths3[$which_path3][2]), 2, 1))
		 {
		 $reconstr = substr($reconstr, 0, 6).substr(($paths3[$which_path3][1]), 2, 1).$state;
		 }

   if (exists $newlist{$paths3[$which_path3][1]}) {}
   else{
   	   $newlist{$paths3[$which_path3][1]} = 0;
	   push @mut, $state;	  
	   #print "$paths3[$which_path3][1] (aa $aa_anc[$i]) to $paths3[$which_path3][2] (aa $aa_des[$i]) is a ", $state, " mut!\n";	  
	   }

# reconstruct the {reconstructed} 1-step codon ($paths3[$which_path3][2]) from the outgroup codon ($codon_OG) 
	   $aa_anc[$i]=codon2aa($codon_OG);
	   $codon_anc[$i]=codonbias($codon_OG);
	   $aa_des[$i]=codon2aa($paths3[$which_path3][2]);
	   $codon_des[$i]=codonbias($paths3[$which_path3][2]);
	   if ($aa_anc[$i] eq $aa_des[$i])
	         {
			 if   (($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'P')) {$state ='SP';}
			 elsif(($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'U')) {$state ='SD';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'P')) {$state ='SB';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'U')) {$state ='SU';} 
			 else  {$state ='S_';}
			}
	   else  {$state = 'R_';}
	   if (substr(($paths3[$which_path3][2]), 0, 1) ne substr($codon_OG, 0, 1))
		 {
		 $reconstr = substr(($paths3[$which_path3][2]), 0, 1).$state.substr($reconstr, 3, 6);
		 }
	   if (substr(($paths3[$which_path3][2]), 1, 1) ne substr($codon_OG, 1, 1))
		 {
		 $reconstr = substr($reconstr, 0, 3).substr(($paths3[$which_path3][2]), 1, 1).$state.substr($reconstr, 6, 3);
		 }
	   if (substr(($paths3[$which_path3][2]), 2, 1) ne substr($codon_OG, 2, 1))
		 {
		 $reconstr = substr($reconstr, 0, 6).substr(($paths3[$which_path3][2]), 2, 1).$state;
		 }

	   
   if (exists $newlist{$paths3[$which_path3][2]}) {}
   else{
	   $newlist{$paths3[$which_path3][2]} = 0;
	   push @mut, $state;	  
	   #print "$paths3[$which_path3][2] (aa $aa_anc[$i]) to $codon_OG (aa $aa_des[$i]) is a ", $state, " mut!\n";	  
	   }
	   
	  #******** push the reconstructed codon to the array called @resconstr_codons
	  push @reconstr_codons, $reconstr;
	#  print "reconstructed codon3: $reconstr\n";

   }
# print "number of mutations only 3 steps= ", (scalar(keys %newlist)-1), "\tcodons: ", (join " ", keys %newlist), "\n";
###########################################################################

# go through all 2-step codons and reconstruct ancestor
for $which_codon2 (0..($num_2-1))   # 2-step mutations 
   {
   $which_path2 = $tables[$i][$which_codon2+$num_3];
   $which_path2 = $which_path2 + $which_codon2*2;
   $reconstr ='NNNNNNNNN';

# reconstruct the 2-step codon ($paths2[$which_path2][0]) from the 1-step codon in its path ($paths2[$which_path2][1]) 
	   $aa_anc[$i]=codon2aa($paths2[$which_path2][1]);
	   $codon_anc[$i]=codonbias($paths2[$which_path2][1]);
	   $aa_des[$i]=codon2aa($paths2[$which_path2][0]);
	   $codon_des[$i]=codonbias($paths2[$which_path2][0]);
	   if ($aa_anc[$i] eq $aa_des[$i])
	         {
			 if   (($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'P')) {$state ='SP';}
			 elsif(($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'U')) {$state ='SD';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'P')) {$state ='SB';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'U')) {$state ='SU';} 
			 else  {$state ='S_';}
			}
	   else  {$state = 'R_';}

	   if (substr(($paths2[$which_path2][0]), 0, 1) ne substr(($paths2[$which_path2][1]), 0, 1))
		 {
		 $reconstr = substr(($paths2[$which_path2][0]), 0, 1).$state.substr($reconstr, 3, 6);
		 }
	   if (substr(($paths2[$which_path2][0]), 1, 1) ne substr(($paths2[$which_path2][1]), 1, 1))
		 {
		 $reconstr = substr($reconstr, 0, 3).substr(($paths2[$which_path2][0]), 1, 1).$state.substr($reconstr, 6, 3);
		 }
	   if (substr(($paths2[$which_path2][0]), 2, 1) ne substr(($paths2[$which_path2][1]), 2, 1))
		 {
		 $reconstr = substr($reconstr, 0, 6).substr(($paths2[$which_path2][0]), 2, 1).$state;
		 }

   if (exists $newlist{$paths2[$which_path2][0]}) {}
   else{
	   $newlist{$paths2[$which_path2][0]} = 0;
	   push @mut, $state;	  
	   #print "$paths2[$which_path2][0] (aa $aa_anc[$i]) to $paths2[$which_path2][1] (aa $aa_des[$i]) is a ", $state, " mut!\n";	  
	   }
	   
# reconstruct the {reconstructed} 1-step codon ($paths2[$which_path2][1]) from the outgroup codon ($codon_OG) 
	   $aa_anc[$i]=codon2aa($codon_OG);
	   $codon_anc[$i]=codonbias($codon_OG);
	   $aa_des[$i]=codon2aa($paths2[$which_path2][1]);
	   $codon_des[$i]=codonbias($paths2[$which_path2][1]);
	   if ($aa_anc[$i] eq $aa_des[$i])
	         {
			 if   (($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'P')) {$state ='SP';}
			 elsif(($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'U')) {$state ='SD';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'P')) {$state ='SB';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'U')) {$state ='SU';} 
			 else  {$state ='S_';}
			}
	   else  {$state = 'R_';}

	   if (substr(($paths2[$which_path2][1]), 0, 1) ne substr($codon_OG, 0, 1))
		 {
		 $reconstr = substr(($paths2[$which_path2][1]), 0, 1).$state.substr($reconstr, 3, 6);
		 }
	   if (substr(($paths2[$which_path2][1]), 1, 1) ne substr($codon_OG, 1, 1))
		 {
		 $reconstr = substr($reconstr, 0, 3).substr(($paths2[$which_path2][1]), 1, 1).$state.substr($reconstr, 6, 3);
		 }
	   if (substr(($paths2[$which_path2][1]), 2, 1) ne substr($codon_OG, 2, 1))
		 {
		 $reconstr = substr($reconstr, 0, 6).substr(($paths2[$which_path2][1]), 2, 1).$state;
		 }

   if (exists $newlist{$paths2[$which_path2][1]}) {}
   else{
   	   $newlist{$paths2[$which_path2][1]} = 0;
	   push @mut, $state;	  
	   #print "$paths2[$which_path2][1] (aa $aa_anc[$i]) to $codon_OG (aa $aa_des[$i]) is a ", $state, " mut!\n";	  
	   }

	  #******** push the reconstructed codon to the array called @resconstr_codons
	  push @reconstr_codons, $reconstr;
	#  print "reconstructed codon2: $reconstr\n";
   }
   
# print " number of mutations incl. 2steps = ", (scalar(keys %newlist)-1), "\tcodons: ", (join " ", keys %newlist), "\n";
###########################################################################

# go through all 1-step codons and reconstruct ancestor
for $which_codon1 (0..($num_1-1)) 
   {
   $reconstr ='NNNNNNNNN';

# reconstruct the 1-step codons ($paths1[$which_codon1]) from the outgroup codon ($codon_OG) 
   $aa_anc[$i]=codon2aa($codon_OG);
   $codon_anc[$i]=codonbias($codon_OG);
   $aa_des[$i]=codon2aa($paths1[$which_codon1]);
   $codon_des[$i]=codonbias($paths1[$which_codon1]);
	   if ($aa_anc[$i] eq $aa_des[$i])
	         {
			 if   (($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'P')) {$state ='SP';}
			 elsif(($codon_anc[$i] eq 'P') and ($codon_des[$i] eq 'U')) {$state ='SD';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'P')) {$state ='SB';}
			 elsif(($codon_anc[$i] eq 'U') and ($codon_des[$i] eq 'U')) {$state ='SU';} 
			 else  {$state ='S_';}
			}
	   else  {$state = 'R_';}

	   if (substr(($paths1[$which_codon1]), 0, 1) ne substr($codon_OG, 0, 1))
		 {
		 $reconstr = substr(($paths1[$which_codon1]), 0, 1).$state.substr($reconstr, 3, 6);
		 }
	   if (substr(($paths1[$which_codon1]), 1, 1) ne substr($codon_OG, 1, 1))
		 {
		 $reconstr = substr($reconstr, 0, 3).substr(($paths1[$which_codon1]), 1, 1).$state.substr($reconstr, 6, 3);
		 }
	   if (substr(($paths1[$which_codon1]), 2, 1) ne substr($codon_OG, 2, 1))
		 {
		 $reconstr = substr($reconstr, 0, 6).substr(($paths1[$which_codon1]), 2, 1).$state;
		 }

   if (exists $newlist{$paths1[$which_codon1]}) {}
   else{
	   $newlist{$paths1[$which_codon1]} = 0;
	   push @mut, $state;	  
	   #print "$paths1[$which_codon1] (aa $aa_anc[$i]) to $codon_OG (aa $aa_des[$i]) is a ", $state, " mut!\n";	  
	   }

	  #******** push the reconstructed codon to the array called @resconstr_codons
	  push @reconstr_codons, $reconstr;
	 # print "reconstructed codon1: $reconstr\n";
   }
# print " number of mutations incl. 1steps = ", (scalar(keys %newlist)-1), "\tcodons: ", (join " ", keys %newlist), "\n";

@path_len[$i] = (scalar(keys %newlist)-1);
# print "length of path $i: $path_len[$i] steps ", @mut, "\n";


# make an array @shortest_path that contains the path (or paths) with the fewest number of steps
if ((scalar(keys %newlist)-1) < $min_length)
   {
   @shortest_path = ();
   $j=0;
   $shortest_path[$j] = (scalar(keys %newlist)-1)." ".(join "", @mut)." ".$i; 
#   print "shortest path was shorter: $shortest_path[$j]\n" ;
   $min_length = (scalar(keys %newlist)-1);
   $countS = ($shortest_path[$j] =~ tr/S//);
   $maxS=$countS;
     }
elsif ((scalar(keys %newlist)-1) == $min_length)
	{
#	print "shortest path was equal:", $shortest_path[$j], "\n" ;
#	print "maxS: $maxS\n";	
	$temp_string = (join "", @mut);
	$countS = ($temp_string =~ tr/S//);
	if ($countS > $maxS)
		{
		$j=0;
		@shortest_path = ();
		$shortest_path[$j] = (scalar(keys %newlist)-1)." ".(join "", @mut)." ".$i; 
		$maxS=$countS;
		}
	elsif ($countS == $maxS)
		{
		$j++;
		$shortest_path[$j] = (scalar(keys %newlist)-1)." ".(join "", @mut)." ".$i;
#		print $shortest_path[$j], "\n" ;
		}
	  }
# print "countS: ",$maxS, 	"\n\n" ;  
}
#######end of the big loop#########

# print the @shortest_path array
# print "shortest paths\n";
#	for $k (0..(scalar @shortest_path-1))
#		{
#		print $shortest_path[$k], "\t" ;
#		print "countS: ", ($shortest_path[$k] =~ tr/S//), "\n";		
#		}

# print reconstructed codons

# some logic: if multiple paths of equal length and equal number of synonymous changes
# we maximise the number of segregating synonymous mutations
# e.g. try AAA GGG TTT ATG
# there are two equally parsimonious paths back to the ancestor that maximise the 
# number of synonymous sites: 
# 1. ATG -> AAG(R) -> AAA(S) -- i.e. the G poly is synonymous here
# 2. ATG -> ATA(R) -> AAA(R) -- i.e. the G poly is a replacement here

$num_codons=$num_1+$num_2+$num_3;
$countS2_max=-1;

# print "reconstructed codons: \n";
for $a (0..(scalar @shortest_path-1))
	{
	@min_paths = ();
	($length, $path, $path_number) = split(/ /, $shortest_path[$a]);
#	print "path number: $path_number\n";
	$countS2=0;
	for $i (($path_number*$num_codons)..((($path_number*$num_codons)+$num_codons)-1))
		{
#		print $reconstr_codons[$i], "\t";
		push @min_paths, $reconstr_codons[$i];
		$countS2 += ($reconstr_codons[($i)] =~ tr/S//);
		}
#	print "\n";
	
	if ($countS2 > $countS2_max)
		{
		@final_codons =();
		@final_codons = @min_paths;
		$log_number = $path_number;
		$countS2_max = $countS2;
		}

	# print combinations for shortest path to check that this is right
#	for $y (0..($num_3+$num_2-1))
#		{
#		print $tables[$path_number][$y] , " ";
#		}
#		print "\n";
	}


# print "best path: ", join(" ", @final_codons), "\tpath number: ", $log_number, "\n";



# make a hash containing each codon, that point to their reconstructed states
# for $num_codons
%final_translations = ();
$final_translations{$codon_OG}='NNNNNNNNN';
for $i (0..($num_codons-1))
	{
	$final_key = substr($final_codons[$i], 0, 1).substr($final_codons[$i], 3, 1).substr($final_codons[$i], 6, 1);
	if (substr($final_codons[$i], 0, 1) eq 'N')
		{
		$final_key = (substr($codon_OG, 0, 1)).(substr($final_key, 1, 2));
		}
	if (substr($final_codons[$i], 3, 1) eq 'N')
		{
		$final_key = (substr($final_key, 0, 1)).(substr($codon_OG, 1, 1)).(substr($final_key, 2, 1));
		}
	if (substr($final_codons[$i], 6, 1) eq 'N')
		{
		$final_key = (substr($final_key, 0, 2)).(substr($codon_OG, 2, 1));
		}
	# print keys to the hash	
	$final_translations{$final_key}=$final_codons[$i];
	
	}

# print "hash final translations\n";
# for $a (keys(%final_translations))
#	{
#	print "key ",  $a, "\t", $final_translations{$a}, "\n";
#	}

push @final_hashes, {%final_translations}; 

} # end loop out


#***** print array of hashes ***************************
# for $i (0..$#final_hashes)
#  {
#   print "$i is { ";
#   for $codon (keys % {$final_hashes[$codon]})
#		{
#        print "$codon=$final_hashes[$i]{$codon} ";
#		}
#		print "}\n";
#  }
   
#********** make strings of paths and minimise R ***********************
# make array of hashes  @final_hash_minR
 @final_hash_minR = ();
 $countR_min = $num_3*3+$num_2*2+$num_1;
# print "max R ", $countR_min, "\n";
 for $i (0..$#final_hashes)
  {
	 @string_array = ();
     for $codon (keys % {$final_hashes[$codon]})
		{
        push @string_array, $final_hashes[$i]{$codon};
		}
	 $string = join ("", @string_array);
	 $countR = ($string =~ tr/R//);
 	 # print "num Rs: $countR\n";
	 if ($countR < $countR_min)
		{
		@final_hash_minR = ();
		push @final_hash_minR, {%{$final_hashes[$i]}};
		$countR_min = $countR;
		}
	 elsif ($countR == $countR_min)
		{
		push @final_hash_minR, {%{$final_hashes[$i]}};
		}
  }
  
if ($complexity eq 'ja')
   {
	print OUT_DIFF "\nHash final_hash_minR no. of elements: ", $#final_hash_minR, "\n";

	for $i (0..$#final_hash_minR)
		{
		print OUT_DIFF "$i is { ";
		for $codon (keys % {$final_hash_minR[$codon]})
			{
			print OUT_DIFF "$codon=$final_hash_minR[$i]{$codon} ";
			}
		print OUT_DIFF "}\n";
		}
	}


#*********************************
# take original codons (@codon_list) and translate them into 6-character code (@codon_list_mod)
#



@codon_list_mod = ();

for $i (0..$#final_hash_minR)
   {
   for ($j=0; $j<(scalar @codon_list_pop); $j++)
      {
      push @codon_list_mod , ($final_hash_minR[$i]{$codon_list_pop[$j]});
#	  $codon_list_mod[$j] = ($final_hash_minR[$i]{$codon_list_pop[$j]});
#	  print  $codon_list_mod[$j], " inside loop\n";
	 }
#	print join ("\-", @{$codon_list_mod_2dim[$i]}), "\n";
#	print join ("-", @codon_list_mod), "\n";

    } 

## return data

 if (scalar(@codon_list_mod)<2)
	{
	print "error in codon_list_mod\n";
	return 0;
    }
else 
    {
#	print join ("\t", @codon_list_mod), "\n";
	return (@codon_list_mod);
	}
 
 
 }



# *************************************************************************
# SUBROUTINE combiner
# 
# input is no. of codon wich are two and three steps away from ancestral codon
# the program returns an array with all possible comparisons for these codons

sub combiner {

($num_3step, $num_2step) = @_;

$num_3comp = (6 ** $num_3step);
$num_2comp = (2 ** $num_2step);
$num_comp = $num_3comp * $num_2comp;

 @comb3=();
 @comb3_dup=();

# first fill out the array with 3 step codons

if ($num_3step>0)
{
for $y (0..($num_3step-1))
  { 
  for $x (0..($num_3comp-1))
	 {
	 $comb3[$x][$y] = int($x%6**($y+1)/(6**$y));
	 }
  }
 }
 else 
 {
 $comb3[0][0]=0;
 $comb3[1][0]=1;
 }
# for each 2step codon we want 
# 1. to duplicate the existing array
# 2. add a column and put "0" in the last column of the first array copy and "1"likewise into the duplicate array
# 3. combine the two arrays
 
for $a (0..($num_2step-1))
	{
	$num_rows = @comb3;
	$num_cols= @{@comb3[0]};

# duplicate the array for later
	@temp = ();
	$x=$y=0;
	for ($x=0; $x<$num_rows; $x++)
		{
		for ($y=0; $y<$num_cols; $y++)
			{
			$temp[$x][$y] = $comb3[$x][$y];
			}
		}
	@comb3_dup = @temp;

	$i=0;
	for ($i=0; $i<$num_rows; $i++)
		{
		$comb3[$i][$num_cols] = 0;         # add "0" to one of them
		$comb3_dup[$i][$num_cols] = 1;    # add "1" to the other
		}

	push (@comb3, @comb3_dup); # combine the two arrays
	}

return (@comb3);
}


#**************************************************************************************************************
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aa {
    my($codon_look) = @_;

    $codon_look = uc $codon_look;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon_look}) {
        return $genetic_code{$codon_look};
    }else{
			return "gap";
            # print STDERR "Bad codon \"$codon_look\"!!\n";
            # exit;
    }
}

       
       
       
  #**************************************************************************************************************
# synNonsyn sites
#
# A subroutine to count the number of synonymous and nonsynonymous sites
# , using hash lookup. The subroutine returns the number of synonymous sites for a codon.
# The number of nonSyn sites = 3 - nr.syn sites

sub countSyn {
    my($codon_look) = @_;

    $codon_look = uc $codon_look;
 
    my(%genetic_code) = (
    
    'TCA' => 1,    # Serine
    'TCC' => 1,    # Serine
    'TCG' => 1,    # Serine
    'TCT' => 1,    # Serine
    'TTC' => (1/3),    # Phenylalanine
    'TTT' => (1/3),    # Phenylalanine
    'TTA' => (2/3),    # Leucine
    'TTG' => (2/3),    # Leucine
    'TAC' => 1,    # Tyrosine
    'TAT' => 1,    # Tyrosine
    #'TAA' => '_',    # Stop
    #'TAG' => '_',    # Stop
    'TGC' => 0.5,    # Cysteine
    'TGT' => 0.5,    # Cysteine
    #'TGA' => '_',    # Stop
    'TGG' => 0,    # Tryptophan
    'CTA' => (4/3),    # Leucine
    'CTC' => 1,    # Leucine
    'CTG' => (4/3),    # Leucine
    'CTT' => 1,    # Leucine
    'CCA' => 1,    # Proline
    'CCC' => 1,    # Proline
    'CCG' => 1,    # Proline
    'CCT' => 1,    # Proline
    'CAC' => (1/3),    # Histidine
    'CAT' => (1/3),    # Histidine
    'CAA' => (1/3),    # Glutamine
    'CAG' => (1/3),    # Glutamine
    'CGA' => (4/3),    # Arginine
    'CGC' => 1,			# Arginine
    'CGG' => (4/3),    # Arginine
    'CGT' => 1,		# Arginine
    'ATA' => (2/3),    # Isoleucine
    'ATC' => (2/3),    # Isoleucine
    'ATT' => (2/3),    # Isoleucine
    'ATG' => 0,    # Methionine
    'ACA' => 1,    # Threonine
    'ACC' => 1,    # Threonine
    'ACG' => 1,    # Threonine
    'ACT' => 1,    # Threonine
    'AAC' => (1/3),    # Asparagine
    'AAT' => (1/3),    # Asparagine
    'AAA' => (1/3),    # Lysine
    'AAG' => (1/3),    # Lysine
    'AGC' => (1/3),    # Serine
    'AGT' => (1/3),    # Serine
    'AGA' => (5/6),    # Arginine
    'AGG' => (2/3),    # Arginine
    'GTA' => 1,    # Valine
    'GTC' => 1,    # Valine
    'GTG' => 1,    # Valine
    'GTT' => 1,    # Valine
    'GCA' => 1,    # Alanine
    'GCC' => 1,    # Alanine
    'GCG' => 1,    # Alanine
    'GCT' => 1,    # Alanine
    'GAC' => (1/3),    # Aspartic Acid
    'GAT' => (1/3),    # Aspartic Acid
    'GAA' => (1/3),    # Glutamic Acid
    'GAG' => (1/3),    # Glutamic Acid
    'GGA' => 1,    # Glycine
    'GGC' => 1,    # Glycine
    'GGG' => 1,    # Glycine
    'GGT' => 1,    # Glycine
    );

    if(exists $genetic_code{$codon_look}) {
        return $genetic_code{$codon_look};
    }else{
			return "gap";
            # print STDERR "Bad codon \"$codon_look\"!!\n";
            # exit;
    }
}

       
       
       
       		  
#**************************************************************************************************************
# codonbias
#
# A subroutine to assing codons as preferred and unpreferred
#   Version 3, using hash lookup

sub codonbias {
    my($codon_bias) = @_;

    $codon_bias = uc $codon_bias;
 
    my(%codon_table) = (
    
    'GCA' => 'U',    # Alanine  ------------------------------------------
    'GCC' => 'P',    # Alanine
    'GCG' => 'P',    # Alanine
    'GCT' => 'U',    # Alanine
    'TGC' => 'P',    # Cysteine  ------------------------------------------
    'TGT' => 'U',    # Cysteine
    'GAC' => 'P',    # Aspartic Acid  ------------------------------------------
    'GAT' => 'U',    # Aspartic Acid
    'GAA' => 'U',    # Glutamic Acid  ------------------------------------------
    'GAG' => 'P',    # Glutamic Acid
    'TTC' => 'P',    # Phenylalanine  ------------------------------------------
    'TTT' => 'U',    # Phenylalanine
    'GGA' => 'U',    # Glycine  ------------------------------------------
    'GGC' => 'P',    # Glycine
    'GGG' => 'P',    # Glycine
    'GGT' => 'U',    # Glycine
	'CAC' => 'P',    # Histidine  ------------------------------------------
    'CAT' => 'U',    # Histidine
    'ATA' => 'U',    # Isoleucine  ------------------------------------------
    'ATC' => 'P',    # Isoleucine
    'ATT' => 'U',    # Isoleucine
    'AAA' => 'U',    # Lysine  ------------------------------------------
    'AAG' => 'P',    # Lysine
    'TTA' => 'U',    # Leucine  ------------------------------------------
    'TTG' => 'U',    # Leucine
    'CTA' => 'U',    # Leucine
    'CTC' => 'P',    # Leucine
    'CTG' => 'P',    # Leucine
    'CTT' => 'U',    # Leucine 
    'AAC' => 'P',    # Asparagine  ------------------------------------------
    'AAT' => 'U',    # Asparagine
    'CCA' => 'U',    # Proline  ------------------------------------------
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'U',    # Proline
    'CAA' => 'U',    # Glutamine  ------------------------------------------
    'CAG' => 'P',    # Glutamine
    'AGA' => 'U',    # Arginine  ------------------------------------------
    'AGG' => 'U',    # Arginine
    'CGA' => 'U',    # Arginine  ------------------------------------------
    'CGC' => 'P',    # Arginine
    'CGG' => 'P',    # Arginine
    'CGT' => 'U',    # Arginine
	'TCA' => 'U',    # Serine   ------------------------------------------
    'TCC' => 'P',    # Serine
    'TCG' => 'P',    # Serine
    'TCT' => 'U',    # Serine
    'AGC' => 'P',    # Serine  ------------------------------------------
    'AGT' => 'U',    # Serine
    'ACA' => 'U',    # Threonine  ------------------------------------------
    'ACC' => 'P',    # Threonine
    'ACG' => 'P',    # Threonine
    'ACT' => 'U',    # Threonine
    'GTA' => 'U',    # Valine  ------------------------------------------
    'GTC' => 'P',    # Valine
    'GTG' => 'P',    # Valine
    'GTT' => 'U',    # Valine
    'TAC' => 'P',    # Tyrosine  ------------------------------------------
    'TAT' => 'U',    # Tyrosine

    'TGA' => '_',    # Stop
    'ATG' => 'M',    # Methionine
    'TGG' => 'W',    # Tryptophan
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop

    );

    if(exists $codon_table{$codon_bias}) {
        return $codon_table{$codon_bias};
    }else{
			return "problem codon bias";
            # print STDERR "Bad codon \"$codon_bias\"!!\n";
            # exit;
    }
}

       
       
   
        
       
       		  
#**************************************************************************************************************
# 4-fold degenerate sites
#
# A subroutine to assing 4-fold degenerate sites
#

sub codon_fourfold {
    my($codon_deg) = @_;

    $codon_deg = uc $codon_deg;
 
    my(%codon_table) = (
    
    'GCA' => 4,    # Alanine  ------------------------------------------
    'GCC' => 4,    # Alanine
    'GCG' => 4,    # Alanine
    'GCT' => 4,    # Alanine
    'TGC' => 2,    # Cysteine  ------------------------------------------
    'TGT' => 2,    # Cysteine
    'GAC' => 2,    # Aspartic Acid  ------------------------------------------
    'GAT' => 2,    # Aspartic Acid
    'GAA' => 2,    # Glutamic Acid  ------------------------------------------
    'GAG' => 2,    # Glutamic Acid
    'TTC' => 2,    # Phenylalanine  ------------------------------------------
    'TTT' => 2,    # Phenylalanine
    'GGA' => 4,    # Glycine  ------------------------------------------
    'GGC' => 4,    # Glycine
    'GGG' => 4,    # Glycine
    'GGT' => 4,    # Glycine
	'CAC' => 2,    # Histidine  ------------------------------------------
    'CAT' => 2,    # Histidine
    'ATA' => 3,    # Isoleucine  ------------------------------------------
    'ATC' => 3,    # Isoleucine
    'ATT' => 3,    # Isoleucine
    'AAA' => 2,    # Lysine  ------------------------------------------
    'AAG' => 2,    # Lysine
    'TTA' => 2,    # Leucine  ------------------------------------------
    'TTG' => 2,    # Leucine
    'CTA' => 4,    # Leucine ------------
    'CTC' => 4,    # Leucine
    'CTG' => 4,    # Leucine
    'CTT' => 4,    # Leucine 
    'AAC' => 2,    # Asparagine  ------------------------------------------
    'AAT' => 2,    # Asparagine
    'CCA' => 4,    # Proline  ------------------------------------------
    'CCC' => 4,    # Proline
    'CCG' => 4,    # Proline
    'CCT' => 4,    # Proline
    'CAA' => 2,    # Glutamine  ------------------------------------------
    'CAG' => 2,    # Glutamine
    'AGA' => 2,    # Arginine  ------------------------------------------
    'AGG' => 2,    # Arginine
    'CGA' => 4,    # Arginine  ------------------------------------------
    'CGC' => 4,    # Arginine
    'CGG' => 4,    # Arginine
    'CGT' => 4,    # Arginine
	'TCA' => 4,    # Serine   ------------------------------------------
    'TCC' => 4,    # Serine
    'TCG' => 4,    # Serine
    'TCT' => 4,    # Serine
    'AGC' => 2,    # Serine  ------------------------------------------
    'AGT' => 2,    # Serine
    'ACA' => 4,    # Threonine  ------------------------------------------
    'ACC' => 4,    # Threonine
    'ACG' => 4,    # Threonine
    'ACT' => 4,    # Threonine
    'GTA' => 4,    # Valine  ------------------------------------------
    'GTC' => 4,    # Valine
    'GTG' => 4,    # Valine
    'GTT' => 4,    # Valine
    'TAC' => 2,    # Tyrosine  ------------------------------------------
    'TAT' => 2,    # Tyrosine

    'TGA' => '_',    # Stop
    'ATG' => 'M',    # Methionine
    'TGG' => 'W',    # Tryptophan
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop

    );

    if(exists $codon_table{$codon_deg}) {
        return $codon_table{$codon_deg};
    }else{
			return "problem with 4-fold degenerate codons";
            # print STDERR "Bad codon \"$codon_deg\"!!\n";
            # exit;
    }
}
       

     
       
   
