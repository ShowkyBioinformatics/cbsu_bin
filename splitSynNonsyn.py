#! /usr/bin/env python

import sys
# import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

### help message
def print_help():
    print "usage: splitSynNonsyn.py <CDS.fasta>"

#read command line return help if only program is called
if len(sys.argv) == 0:
    sys.exit(print_help())
else:
    input_CDS = sys.argv[1]

### input fasta file is only command line argument
# input_CDS = sys.argv[1]

### Specify two output fasta file names
output_Syn = input_CDS.split(".")[0] + ".Syn.fasta"
output_nonSyn = input_CDS.split(".")[0] + ".nonSyn.fasta"

### read in fasta file
seq_records = SeqIO.parse(input_CDS, "fasta")


### nonsynonymous sites extraction function
def getNonsyonymousSites(records):
    """ Extract 1st and 2nd nucleotide for each 
        codon in CDS sequence"""
    for record in records:
        nonSyn = Seq("")
        for i in range(0, len(record), 3):
            nonSyn += record[i:i+2:]
        yield nonSyn

### extract nonsynonymmous sites
nonSynSites = getNonsyonymousSites(seq_records)

### Output nonsynonymous sites fasta file and print message
count = SeqIO.write(nonSynSites, output_nonSyn,"fasta")
print("extracted nonsynonymous sites from %i sequences" %count)

### read in fasta file
seq_records = SeqIO.parse(input_CDS, "fasta")

### synonymous sites extraction function
def getSynonymousSites(records):
    """ Extract 3rd nucleotide for each codon in 
        CDS sequence"""
    for record in records:
        yield record[2::3]

### extract synonymous sites
synSites = getSynonymousSites(seq_records)

### Output nonsynonymous sites fasta file and print message
count = SeqIO.write(synSites,output_Syn,"fasta")
print("extracted synonymous sites from %i sequences" %count)
