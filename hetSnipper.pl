#!/usr/bin/perl

=head1 SYNOPSIS

extract heterospecific SNPs:

hetSnipper.pl <gene_list> <BAM_list> <reference> <output_folder>

  
=cut

# required perl modules
use strict;
use warnings;
use Pod::Usage;

# assign command line arguments to variables
my $gene_list = $ARGV[0];
my $ref = $ARGV[2];
my $output = $ARGV[3];
my ($help);

pod2usage(1) if($help);
pod2usage("No files given!\n")  if ((!$gene_list));

# read gene list file
open (FILE1, $gene_list) or die ("Could not open gene list file \n");

# create output directory
system ("mkdir $output");

# start a while loop with the gene list
while ($gene_list = <FILE1>){

	# read genes name one by one and assign them to an array
	chomp $gene_list;
	my @gene_bits = split(/\s+/,$gene_list);
	my $gdir_name = $gene_bits[0];
	my @garray = ();
	push (@garray, $gdir_name);

	# loop over each gene in the array
	foreach my $gint (@garray){

		# assign the bam file from the command line to a variable
		my $BAM_list = $ARGV[1];

		# read the bam list file
		open (FILE2, $BAM_list) or die ("Could not open gene list file \n");

		# create output directory for current gene
		system ("mkdir $output/$gint.output");

		# output fasta file for current gene
		system ("faOneRecord $ref $gint > $gint.fasta");
		system ("mv $gint.fasta $output/$gint.output");

		# start a while loop with the bam list
		while($BAM_list = <FILE2>){

			# read bam file names one by one and assign them to an array
			chomp $BAM_list;
			my @bam_bits = split(/\s+/,$BAM_list);
			my $bdir_name = $bam_bits[0];
			my @barray = ();
			push (@barray, $bdir_name);

			# loop over each bam file in the array
			foreach my $bint (@barray){

				# output bam file for the current gene
				system ("samtools view -b $bint $gint > $gint.$bint.bam");

				# output a VCF file for current gene alignment
				system ("samtools mpileup -o $gint.$bint.vcf -v -u -f $ref -r $gint $bint");

				# filter out SNPs by depth (here by 50 reads)
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
