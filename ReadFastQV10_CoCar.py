import os
import gzip
import time
import glob
import argparse
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import warnings
warnings.filterwarnings('ignore')



def get_arguments():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('-R', '--runName', dest='runName', type=str,
                        help='Name of run')
    parser.add_argument('-O', '--output', dest='output', type=str,
                        help='/path/to/output')
    parser.add_argument('-I', '--input', dest='input', type=str,
                        help='/path/to/fastq files')
    parser.add_argument('-T', '--threads', dest='threads', type=str, default=os.cpu_count(),
                        help='number of threads to use for paralleling (default maximum cpu threads)')
    return parser.parse_args() 


def read_fastq(fastq, fastq_dir, out_directory):
    step = 4
    reads = []
    n = 0
    #print('Opening_FastQ')
    #print(fastq)
    os.chdir(fastq_dir)
    with gzip.open(fastq) as file:
        for lineno, line in enumerate(file):
            if (lineno-1) % step == 0:
                reads.append(line)
                n += 1
                #if n % 1000000 == 0:
                    #print('  ' + str(round(n / 1000000)) + ' M Reads')
    #print(str(len(reads)) + '  Reads')
    #print(' ')
    os.chdir(out_directory)
    filename = str(fastq.split("_")[0]) + '_Reads.csv'
    np.savetxt(filename, reads, delimiter=";", fmt='%s')



def main():
    start_time = time.time()

    #get arguments
    args = get_arguments()
    out_directory = os.path.realpath(args.output) + os.path.sep + args.runName
    fastq_dir = os.path.realpath(args.input)
    threads_number = args.threads
    
    #create directory if it not exists
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    #get sample list
    os.chdir(fastq_dir)
    fastq_list = glob.glob("*.fastq.gz")

    #launch analysis
    Parallel(n_jobs= threads_number, verbose = 1)(delayed(read_fastq)(fastq, fastq_dir, out_directory) for fastq in fastq_list)

    print(f"--- Program executed in {float(time.time() - start_time) / 60} minutes ---")



if __name__ == "__main__":
    main()
