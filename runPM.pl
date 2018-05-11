#!/usr/bin/env perl

########################
### Required Modules ###

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use FindBin;
use File::Basename;
use Bio::Tools::Run::Alignment::Clustalw;
# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);
# for input of the sequence data
use Bio::SeqIO;
use Bio::AlignIO;

my $usage = <<__EOUSAGE__;

#######################################################################################################################################################################
#   
#    ____       _                                  _                           
#   |  _ \\ ___ | |_   _ _ __ ___   ___  _ __ _ __ | |__   ___  _ __ ___   __ _ 
#   | |_) / _ \\| | | | | '_ ` _ \\ / _ \\| '__| '_ \\| '_ \\ / _ \\| '_ ` _ \\ / _` |
#   |  __/ (_) | | |_| | | | | | | (_) | |  | |_) | | | | (_) | | | | | | (_| |
#   |_|   \\___/|_|\\__, |_| |_| |_|\\___/|_|  | .__/|_| |_|\\___/|_| |_| |_|\\__,_|
#                 |___/                     |_|                                
#
#   (Poymorphoma was written by Doris Bachtrog and Peter Andolfatto)
#
#
#       Requirements:
#           
#           1. In addition to the perl modules required for this script to work (see code header), the following scripts are required:
#               
#               - Perl scripts:
#
#                   fastagrep.pl
#                   MK.pl   
#                   MK_trim.pl  
#
#           2. CDS file should contain all sequences to be analyzed. Transcript ID within species (among lines) should be the same. Header should have follwoing format:
#               >transcript_id .......... species=<species_id> line=<line_id>
#
#           
#   Required arguments:
#
#   --CDS_file|c <string>               sequence file containing CDS from multiple species/lines.
#
#   --gene_list|g <string>              list file of genes to be analyzed.
#       OR
#   --gene_id|i <string>                transcript/gene to be analyzed
#   
#   --output|o                  name of directory to place outputs 
#
#
#
#   Options:
#
#   --show_samples                  display species/line IDs in CDS_file and exit (only requires "-c <CDS_file>" argument).
#
#   --include_samples               include specified species/lines 
#
#       --samples_file <string>         (use line ID if present, otherwise use species ID).
#
#   --outgroup_id             
#
#   --view_DNA_alignment                open DNA alignment in Geneious.
#
#   --view_Protein_alignment            open Protein alignment in Geneious.
#
#   --save_alignments               retain alignment files
#
#
#   --PM_coding
#
#   --PM_noncoding
#
#
#   --help                      
#
#                           
#######################################################################################################################################################################


__EOUSAGE__

    ;


#####################################
## Set BioPerl alignment protocol ###

my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new;


###################################
### Define command line options ###

my $CDS_file;
my $gene_list;
my $gene_id;

my $output_dir;

my $show_samples_flag = 0;
my $include_samples_flag = 0;
my $samples_file;

my $outgroup_id;

my $view_aln_dna_geneious_flag = 0;
my $view_aln_aa_geneious_flag = 0;
my $save_alignments_flag = 0;

my $PM_coding = 0;
my $PM_noncoding = 0;

my $help_flag = 0;


&GetOptions (  

    'CDS_file|c=s' => \$CDS_file,
    'gene_list|g=s' => \$gene_list,
    'gene_id|i=s' => \$gene_id,

    'output|o=s' => \$output_dir,
    
    'show_samples' => \$show_samples_flag,
    'include_samples' => \$include_samples_flag,
    'samples_file=s' => \$samples_file,
    
    'outgroup_id=s' => \$outgroup_id,


    'view_DNA_alignment' => \$view_aln_dna_geneious_flag,
    'view_Protein_alignment' => \$view_aln_aa_geneious_flag,
    'save_alignments' => \$save_alignments_flag,
    
    'PM_coding' => \$PM_coding,
    'PM_noncoding' => \$PM_coding,
    
    'help' => \$help_flag,
    
);


#####################################################
### Check command line arguments and housekeeping ###

# if command line argument are not recognized
if (@ARGV) {
    die "Error, don't understand parameters: @ARGV";
}

if ($help_flag) {
    die $usage;
}

# only show samples in CDS_file
if ($show_samples_flag) {
    system ("grep '>' $CDS_file | sed 's/.*species=/species=/g' | sort -u");
    exit(0);
}

# if Mandatory arguments not present
unless ($gene_list || $gene_id || $output_dir) {
    die "\nNo CDS file, gene ID or gene list file provided.\n\nAnd you need to specify output directory!\n $usage";
}

if (@ARGV) {
    die "Error, do not recognize params: @ARGV ";
}

# specify output dir cname
if ($output_dir) {
    mkdir($output_dir);
}

# specify output dir cname
if ($PM_coding || $PM_noncoding) {
    system ("mkdir nex.files");
}

# if Mandatory arguments not present
if ($PM_coding || $PM_noncoding) {
        unless ($outgroup_id){
            die "\nNo outgroup specified.\n\nNeed an outgroup for Polymorphoma analyses\n $usage";
        }
}



######################################
###### Main execution of processes ###

main: {
                
        my @gene_object;
        
        ###############################################################
        ##### Set the gene object to a single gene or list of genes. ##
        if ($gene_id) {
            push (@gene_object, $gene_id);
        } elsif ($gene_list) {
            open (FILE1, $gene_list) or die ("Could not open file \n");
            while ($gene_list = <FILE1>){
                chomp $gene_list;
                my @file_bits = split(/\s+/,$gene_list);
                my $gene_names = $file_bits[0];
                push (@gene_object, $gene_names);
            }
        }
        
        
        ###############################################################
        ### Initiate the loop to analyze individual transcripts #######
    
        foreach my $transcript (@gene_object){
            eval {
            
            ###########################################################
            ### Add orthologous sequences to the analysis...  #########
            ### Processes FASTA header to include species and #########
            ### line ID, if present.                          #########
            
            system ("fastagrep.pl -X $transcript $CDS_file > $transcript.all.fa");

            if ($include_samples_flag) {
                system ("fastagrep.pl -f $samples_file $transcript.all.fa > $transcript.sample_subset.fa");
                system ("rm $transcript.all.fa");
                system ("mv $transcript.sample_subset.fa $transcript.fa");
                # header processing, contd.
                system ("sed -i.bak '/>/ s/ line=/_/g' $transcript.fa");
                system ("sed -i.bak '/>/ s/>.*species=/>/g' $transcript.fa");
                }
            else {
                system ("sed -i.bak '/>/ s/ line=/_/g' $transcript.all.fa");
                system ("sed '/>/ s/>.*species=/>/' $transcript.all.fa > $transcript.fa");
                system ("rm $transcript.all.fa");
            }

            # read in sequence  
            my $seqio = Bio::SeqIO->new(-file => "$transcript.fa",
                                       -format => 'fasta');
            my %seqs;
            my @prots;
            
            # process each sequence
            while ( my $seq = $seqio->next_seq ) {
                $seqs{$seq->display_id} = $seq;
                
                # translate them into protein
                my $protein = $seq->translate();
                my $pseq = $protein->seq();
                if( $pseq =~ /\*/ &&
                $pseq !~ /\*$/ ) {
                    warn("Oi! $transcript has a stop codon!!!");
                }
                
                # Tcoffee can't handle '*' even if it is trailing
                $pseq =~ s/\*//g;
                $protein->seq($pseq);
                push @prots, $protein;
            }
            
            # Warn if only 1 sequence present
            if( @prots < 2 ) {
                warn("$transcript - Need at least 2 CDS sequences to proceed");
            }

            # Align the sequences with clustalw
            my $aa_aln = $aln_factory->align(\@prots);
            
            # project the protein alignment back to CDS coordinates
            my $dna_aln = aa_to_dna_aln($aa_aln, \%seqs); 
            my @each = $dna_aln->each_seq();
            
            # output alignments for downstream analysis
            my $out_dna = Bio::AlignIO->new(-file => ">$transcript.aln_dna.afa" ,
                                        -format => 'fasta');
            $out_dna -> write_aln($dna_aln);

            my $out_aa = Bio::AlignIO->new(-file => ">$transcript.aln_aa.afa" ,
                                        -format => 'fasta');
            $out_aa -> write_aln($aa_aln);
            
            # clean-up alignment header files
            system ("sed -i.bak '/>/ s/\\/.*//g' $transcript.*.afa");
            system ("rm *bak");
            
            
            ############################################################
            ### Perform the analyses specified in the command line #####
            
            # MK test
            if ($PM_coding) {
            
                # Edit/move alignment file to MK tmp dir
                system ("faOneRecord $transcript.aln_dna.afa $outgroup_id > $transcript.aln_dna.ord.afa");
                system ("fastagrep.pl -X -v $outgroup_id $transcript.aln_dna.afa >> $transcript.aln_dna.ord.afa");
                system ("mv $transcript.aln_dna.ord.afa $transcript.nex");
                system ("mv $transcript.nex nex.files");
                # system ("rm $transcript.tmp.nex");
                
                
                # Clean-up seqs from MK tmp directory for the next sequence
                # system ("rm MK.tmp/$transcript.fa ");
            }   
            
            
                
            
            
            # option to save the alignment files
            if ($save_alignments_flag || $view_aln_dna_geneious_flag || $view_aln_aa_geneious_flag) {
                system ("mv $transcript.aln_dna.afa $output_dir/$transcript.DNA_alignment.afa");
                system ("mv $transcript.aln_aa.afa $output_dir/$transcript.Protein_alignment.afa");
            } else {
                system ("rm $transcript.aln*.afa");
            }
            
            # option to view alignment files in Geneious
            if ($view_aln_dna_geneious_flag) {
                system ("open -a Geneious $output_dir/$transcript.DNA_alignment.afa");
            }
                
            
            if ($view_aln_aa_geneious_flag) {
                system ("open -a Geneious $output_dir/$transcript.Protein_alignment.afa");
            }
            
            
            # clean-up
            system ("rm $transcript.fa");
            
            ##########################################################
            ############### END OF LOOP ##############################
        };       
    }

       if ($PM_coding) {
            system ("ls -1 nex.files/* > transcript.list");
            system ("Polymorphorama_coding.pl transcript.list 0.1");
            system ("rm -r nex.files");
            system ("mv summarystats.xls frequencies.xls codon_bias.xls $output_dir");
            system ("rm transcript.list")
        } 
            
}



