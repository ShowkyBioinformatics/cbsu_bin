#!/usr/bin/env python3
'''
Make a samples.json file with sample names and file names (single-end reads).
'''
def msg(name=None):                                                            
    return ''' make_json_SE_samples.py <samples_files>
        '''

import json
from glob import glob
from sys import argv
import argparse
parser = argparse.ArgumentParser(description='Make a samples.json file with sample names and file names (single-end reads).', usage=msg())

# Change this line to match your filenames.
fastqs = glob(argv[1])
FILES = {}

# Change this line to extract a sample name from each filename.
SAMPLES = [fastq.split('/')[-1].split('.')[0] for fastq in fastqs]

for sample in SAMPLES:
    # Change 'R1' to match the way your reads are marked.
    reads = lambda fastq: sample in fastq and 'R1' in fastq
    FILES[sample] = {}
    FILES[sample]['R1'] = sorted(filter(reads, fastqs))

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)

