import os
import shutil
import gzip
import csv
import numpy as np
import pandas as pd

### Procedure generale : recupere les reads d'un FastQ (independant du type d'analyse)

### Nom du run
#run = '2021_04_14_CoCar1'
#run = '2021_06_04_CoCar2'
#run = '2021_09_15_CoCar3'
#run = '2021_09_27_CoCar4'
run = 'Cocar_Test'

rep = 'C:\CoCar\Reads_CoCar'+ os.path.sep + run
os.mkdir(rep)

### Listes FastQ et Echantillons
ListFastQ = []
repFastQ = 'C:\CoCar\FastQ' + os.path.sep + run
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