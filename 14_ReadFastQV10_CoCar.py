import os
import shutil
import gzip
import csv
import numpy as np
import pandas as pd
import argparse

my_parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

my_parser.add_argument('-R','--runName',type=str,
                       help='Name of run')

my_parser.add_argument('-O','--output',type=str,
                       help='/path/to/output')

my_parser.add_argument('-I','--input',type=str,
                       help='/path/to/fastq files')

args = my_parser.parse_args()

run = args.runName

rep = os.path.realpath(args.output)

### Listes FastQ et Echantillons
ListFastQ = []
repFastQ = os.path.realpath(args.input)

rep = rep+ os.path.sep + run
if not os.path.exists(rep):
    os.mkdir(rep)

### Listes FastQ et Echantillons
ListFastQ = []
for file in os.listdir(repFastQ):
    ListFastQ.append(file)

### Script
# recherche des groupes de (4) FastQ
ListSamples = []
for FastQ in ListFastQ:
    #print(FastQ)
    #print(FastQ.find('L00'))
    spl = FastQ[0:FastQ.find('L00')+3]
    if len(ListSamples) == 0 or spl not in ListSamples:
        ListSamples.append(spl)

# Parcours des FastQ - Recuperation des reads
step = 4

for Sample in ListSamples:
    n = 0
    reads = []
    for i in ['1','2','3','4']:
        FastQ = Sample + i + '_R1_001.fastq.gz'
        n = 0
        print('Opening_FastQ')
        print(FastQ)
        os.chdir(repFastQ)
        with gzip.open(FastQ) as f:
            for lineno, line in enumerate(f):
                if (lineno-1) % step ==0:
                    reads.append(line)
                    n=n+1
                    if n % 1000000 == 0:
                        print('  ' + str(round(n / 1000000)) + ' M Reads')
        print(str(len(reads)) + '  Reads')
        print(' ')

    # Creation du repertoire de sauvegarde
    os.chdir(rep)
    namfile = str(Sample[0:4]) + '_Reads.csv'
    np.savetxt(namfile, reads, delimiter=";", fmt='%s')
