#!/usr/bin/env python3
'''
Make a samples.json file with sample names and file names.
'''
def msg(name=None):                                                            
    return ''' make_json_samples.py <samples_files>
        '''

import json
from glob import glob
from sys import argv
import argparse
parser = argparse.ArgumentParser(description='Make a samples.json file with sample names and file names.', usage=msg())

# Change this line to match your filenames.
fastqs = glob(argv[1])
FILES = {}

# Change this line to extract a sample name from each filename.
SAMPLES = [fastq.split('/')[-1].split('.')[0] for fastq in fastqs]

for sample in SAMPLES:
    mate1 = lambda fastq: sample in fastq and 'R1' in fastq
    mate2 = lambda fastq: sample in fastq and 'R2' in fastq
    if any('R2' in s for s in sorted(filter(mate2, fastqs))):
        FILES[sample] = {}
        FILES[sample]['R1'] = sorted(filter(mate1, fastqs))
        FILES[sample]['R2'] = sorted(filter(mate2, fastqs))
    else:
        FILES[sample] = {}
        FILES[sample]['R1'] = sorted(filter(mate1, fastqs))

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)