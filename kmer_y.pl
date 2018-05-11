#!/usr/bin/perl
#kmer.pl
# System Requirements 
# 1. A computer with a Unix/Linux operating system (Tested on CentOS, and Ubuntu) 
# 2. Perl (http://www.perl.org) 
# 3. Bowtie installed and in your path (http://bowtie-bio.sourceforge.net/index.shtml)
#
# To run kmer.pl you need the following data
#1. female kmers produced by jellfish, fasta file
#2. male  kmers produced by jellfish, fasta file
#3. assembled contig sequences (including genome transcript or small RNA sequences), fasta file 
#
#USAGE:
#perl kmer.pl [options]
#Options:
# -k [default parameter. Producing W-specific kmers fasta file. Please use -nok if W-specific kmers fasta file already exist]
# -f [female kmer data produced by jellyfish, fasta file. default -k]
# -m [male kmer data produced by jellyfish, fasta file. default -k]
# -a [minimum number of kmer frequency in female reads, default: 15]
# -b [using bowtie to produce AWK value for each already assembled contig]
# -c [contig file name, fasta file. It can be genome contig sequences, transcript sequences or small RNA sequences]
# -t [thread of using bowtie, default: 8]
# -p [output prefix, default: test]
# -h [print help messeage]
#
#Example:
# perl kmer.pl -f female_counts_dump.fa -m male_counts_dump.fa -c assembled_contig.fa -b -t 8 -a 15 -p silkworm_W 
# File containing results with be silkworm_W  

use Getopt::Long;
my $bowtie="",$kmer="";
my $male,$female,$bowtie,$contigs;
my $alignment_threshold = 15,$prefix="test",$thread=8;
GetOptions( 
    "h"=>\$help,
    "k!" => \$kmer,
	"b"	=> \$bowtie,
	"f=s" => \$female,
	"m=s" => \$male,
    "c:s" => \$contigs,
    "t:i" => \$thread, 
    "a=i" => \$alignment_threshold,
    "p=s" => \$prefix
	);
&help();
&get_w_kmer($male,$female);
&run_bowtie_index();
&run_bowtie();
&awk();
sub help{
    $help_message ="\nUSAGE: perl kmer.pl [options]\n\nOptions:\n-k [default parameter. Producing W-specific kmers fasta file. Please use -nok if W-specific kmers fasta file already exist]\n-f [female kmer data produced by jellyfish, fasta file. default -k]\n-m [male kmer data produced by jellyfish, fasta file. default -k]\n-a [minimum number of kmer frequency in female reads, default: 15]\n-b [using bowtie to produce AWK value for each already assembled contig]\n-c [contig file name, fasta file. It can be genome contig sequences, transcript sequences or small RNA sequences]\n-t [thread of using bowtie, default: 8]\n-p [output prefix, default: test]\n-h [print help messeage]\n\nExample:\nperl kmer.pl -f female_counts_dump.fa -m male_counts_dump.fa -c assembled_contig.fa -b -t 8 -a 15 -p silkworm_W\nFile containing results with be silkworm_W\n";
    if($help){
        print "$help_message\n";
        exit 0;
    }
}

sub get_w_kmer {
    if($kmer)
    {
$which_data = $_[0];
open (MALE,$which_data) or die "Can not open Male_kmer$!\n";
$/=">";
%male_kmer;
while(<MALE>){
        chomp;
       ($fre,$seq)=split(/\n+/,$_,2);
        $male_kmer{$seq}=1;
     
}
$prefix_name=$prefix."_W_specific_kmer.fa";
open(FEMALE,$female) or die "Can not open Female_kmer$!\n";
open (OUT,">$prefix_name") or die "Can not open OUTPUT_FILE$!\n";
$/=">";
my $name=1;     
while(<FEMALE>){
        chomp;
                ($frequence,$id)=split(/\n+/,$_,2);
                if((!$male_kmer{$id})&&($frequence>=$alignment_threshold)){
                        print  OUT ">$name\t$frequence\n$id";
                $name++;}

}
}
else{
}
}
sub run_bowtie_index{
    if($bowtie){
    #make the index with bowtie-build 
    
    print "Running bowtie-build...\n";
    system("mkdir index");
    system("bowtie-build $contigs ./index/contigs");
}
}


sub run_bowtie{
    $prefix_name=$prefix."_W_specific_kmer.fa";
    if(($bowtie)&&(-e $prefix_name)){
    print "Aligning reference to database...\n";
    system ("bowtie -f -a -p $thread -v 0 ./index/contigs  $prefix_name --suppress 1,2,4,5,6,7,8,9>bowtie_result 2> bowtie_log"); 
    system ("sort bowtie_result|uniq -c >bowtie_result2")
    }
    if(! -e $prefix_name){
        print "No W specific kmer file exist, Please check your files! (valid file: *_W_specific_kmer.fa, * means -p parameter)\n";
    }
   }


sub awk{
    if($bowtie){
        open(BOWTIE,"bowtie_result2") or die "Can not open Bowtie_result$!\n";
        %hash;
        while(<BOWTIE>){
            chomp;
            ($start,$num,$id)=split(/\s+/,$_);
            $hash{$id}=$num;
        }
        open(CONTIG,$contigs) or die "Can not open Contigs$!\n";
        $prefix_name2=$prefix."_awk.txt";
        open (OUT2,">$prefix_name2") or die "Can not open OUTPUT_FILE2$!\n";
        print OUT2  "ID\tamount of W kmer\tAWK value\n";
        $/=">";
        while(<CONTIG>){
            chomp;
            ($name,$seq,$end)=split(/\n+/,$_);
            $len=length($seq);
	    if($len!=0){
            $awk=sprintf "%.2f",1000*int($hash{$name})/int($len);
}
            if($hash{$name}){
                print OUT2 "$name\t$hash{$name}\t$awk";            
}

        }
    }

}














