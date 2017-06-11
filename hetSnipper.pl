#!/usr/bin/perl

=head1 SYNOPSIS

extract heterospecific SNPs:

hetSnipper.pl <gene_list> <BAM_list> <reference> <output_folder>

  
=cut

use strict;
use warnings;
use Pod::Usage;

my $gene_list = $ARGV[0];
my $ref = $ARGV[2];
my $output = $ARGV[3];
my ($help);

pod2usage(1) if($help);
pod2usage("No files given!\n")  if ((!$gene_list));

open (FILE1, $gene_list) or die ("Could not open gene list file \n");

system ("mkdir $output");

while ($gene_list = <FILE1>){
	chomp $gene_list;
	my @gene_bits = split(/\s+/,$gene_list);
	my $gdir_name = $gene_bits[0];
	my @garray = ();
	push (@garray, $gdir_name);
	foreach my $gint (@garray){
		my $BAM_list = $ARGV[1];
		open (FILE2, $BAM_list) or die ("Could not open gene list file \n");
		system ("mkdir $output/$gint.output");
		system ("faOneRecord $ref $gint > $gint.fasta");
		system ("mv $gint.fasta $output/$gint.output");
		while($BAM_list = <FILE2>){
			chomp $BAM_list;
			my @bam_bits = split(/\s+/,$BAM_list);
			my $bdir_name = $bam_bits[0];
			my @barray = ();
			push (@barray, $bdir_name);
			foreach my $bint (@barray){
				system ("samtools view -b $bint $gint > $gint.$bint.bam");
				system ("samtools mpileup -o $gint.$bint.vcf -v -u -f $ref -r $gint $bint");
				system ("/programs/vcflib/bin/vcffilter -f \"DP > 50\" $gint.$bint.vcf > $gint.$bint.filtered.vcf");
				system ("awk 'BEGIN{FS=OFS=\"\\t\"}{if (\$5 != \"<*>\") print}' $gint.$bint.filtered.vcf | sed '/\#\#/d' | sed 's/BQB=.*;DP/DP/g' | sed 's/DP=.*QS=/QS=/g' | sed 's/;.*//g' | sed 's/QS=/   /g' | awk '{ gsub(\",\", \" \", \$8) ; print }' | awk '{if (\$8 < 0.9) print \$0}' | tr ' ' '\t' > $gint.$bint.result.vcf");
				system ("mv $gint.$bint.bam $gint.$bint.result.vcf $output/$gint.output");
				system ("rm $gint.$bint.vcf $gint.$bint.filtered.vcf");
			}
		}
	}
}
close FILE1;
close FILE2;
