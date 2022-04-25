import os
import glob
import time
import shutil
import joblib
import subprocess
import argparse
import contextlib
from tqdm import tqdm
from joblib import Parallel, delayed

import ReadFastQV10_CoCar
import Count_CoCar
import Figures_CoCar
import Recherche_Anos_CoCar


def get_arguments():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('-R', '--run_modes', dest='mode', type=str, default="align,count,figure,stat,model", 
                        help='run mode option comma separated, default = align,count,figure,stat,model')
    parser.add_argument('-I', '--input', dest='input', type=str,
                        help='/path/to/fastq files')
    parser.add_argument('-O', '--output', dest='output', type=str,
                        help='/path/to/output')
    parser.add_argument('-P', '--probes', dest='probes', type=str, default="RT-MLPA_scripts/Probes_Design", 
                        help='/path/to/probes ref')
    parser.add_argument('-G', '--genes_file', dest='gene_file', type=str,
                        help='/path/to/genes file')
    parser.add_argument('-T', '--threads', dest='threads', type=str, default=os.cpu_count(),
                        help='number of threads to use for paralleling (default maximum cpu threads)')
    return parser.parse_args() 


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()



def main():
    start_time = time.time()

    print('''\
                                                                                   
         /$$$$$$$  /$$$$$$$$         /$$      /$$ /$$       /$$$$$$$   /$$$$$$ 
        | $$__  $$|__  $$__/        | $$$    /$$$| $$      | $$__  $$ /$$__  $$
        | $$  \ $$   | $$           | $$$$  /$$$$| $$      | $$  \ $$| $$  \ $$
        | $$$$$$$/   | $$    /$$$$$$| $$ $$/$$ $$| $$  RL  | $$$$$$$/| $$$$$$$$
        | $$__  $$   | $$   |______/| $$  $$$| $$| $$      | $$____/ | $$__  $$
        | $$  \ $$   | $$           | $$\  $ | $$| $$      | $$      | $$  | $$
        | $$  | $$   | $$    Louise | $$ \/  | $$| $$$$$$$$| $$      | $$  | $$
        |__/  |__/   |__/           |__/     |__/|________/|__/      |__/  |__/''')


    # get arguments
    print("\n--------CONFIG--------")
    args = get_arguments()
    out_directory = os.path.realpath(args.output) + os.path.sep; print(f"Output directory: {out_directory}") 
    fastq_dir = os.path.realpath(args.input); print(f"FASTQ directory: {fastq_dir}")
    gene_file = os.path.realpath(args.gene_file)
    probes_dir = os.path.realpath(args.probes); print(f"Probes directory: {probes_dir}")
    threads_number = args.threads; print(f"number of threads use: {threads_number}")
    genes_list = []
    with open(gene_file, "r") as file:
        for line in file:
            genes_list.append(line.strip())
    print(f"Gene list: {' '.join(genes_list)}")
    run_mode = args.mode; print(f"Runmode: {run_mode}")
    
    align=False
    count=False
    figure=False
    stat=False
    model=False

    for mode in run_mode.split(","):
        if mode == "align":
            align=True
        if mode == "count":
            count=True
        if mode == "figure":
            figure=True
        if mode == "stat":
            stat=True
        if mode == "model":
            model=True
    

    print("\n--------PROCESSING--------")
    if align:
        print("Step 1: read FASTQ...")
        
        if not os.path.exists(out_directory + "/mySequence"):
            os.makedirs(out_directory + "/mySequence")
        os.chdir(fastq_dir)
        fastq_list = glob.glob("*.fastq.gz")
        with tqdm_joblib(tqdm(desc="Read FASTQ", total=len(fastq_list), ncols=75)) as progress_bar:
            Parallel(n_jobs= threads_number, verbose = 0)(delayed(ReadFastQV10_CoCar.read_fastq)(fastq, fastq_dir, out_directory + "/mySequence") for fastq in fastq_list)
    
    if count:
        print("Step 2: count...")
        if not os.path.exists(out_directory + "/myCount"):
            os.mkdir(out_directory + "/myCount")
        os.chdir(out_directory + "/mySequence")
        reads_list = glob.glob("*.csv")
        samples_list = [read.split("_")[0] for read in reads_list]
        with tqdm_joblib(tqdm(desc="Count", total=len(samples_list), ncols=75)) as progress_bar:
            Parallel(n_jobs= threads_number, verbose = 0)(delayed(Count_CoCar.count_analysis)(sample, genes_list, out_directory + "/myCount", out_directory + "/mySequence", probes_dir) for sample in samples_list)
    
    if figure:
        print("Step 3: make figures...")
        samples_list = [file for file in os.listdir(out_directory + "/myCount")]
        with tqdm_joblib(tqdm(desc="Make figure", total=len(samples_list), ncols=75)) as progress_bar:
            Parallel(n_jobs= threads_number, verbose = 0)(delayed(Figures_CoCar.draw_figure)(sample, genes_list, out_directory + "/myCount", probes_dir) for sample in samples_list)
    
    if stat:
        print("Step 4: research ANOS...")
        if not os.path.exists(out_directory + "/myAnalysis" + os.path.sep + '#Moy'):
            os.makedirs(out_directory + "/myAnalysis" + os.path.sep + '#Moy')
        for gene in genes_list:
            if not os.path.exists(out_directory + "/myAnalysis" + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Counts'):
                os.makedirs(out_directory + "/myAnalysis" + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Counts')
            if not os.path.exists(out_directory + "/myAnalysis" + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices'):
                os.makedirs(out_directory + "/myAnalysis" + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices')
        index_matrix = []
        for file in os.listdir(out_directory + "/myCount"):
            if file != 'Matrice' and file != 'Ano_' and file != '#Moy':
                index_matrix.append(file)
        with tqdm_joblib(tqdm(desc="Research ANOS", total=len(genes_list), ncols=75)) as progress_bar:
            Parallel(n_jobs= threads_number, verbose = 0)(delayed(Recherche_Anos_CoCar.search_anos)(gene, index_matrix, out_directory + "/myAnalysis", "myAnalysis", out_directory + "/myCount", probes_dir) for gene in genes_list)
        
        #copy #Moy folder in myCount
        shutil.copytree(out_directory + "/myAnalysis" + os.path.sep + '#Moy', out_directory + "/myCount" + os.path.sep + '#Moy')
        samples_list = ["#Moy"]
        with tqdm_joblib(tqdm(desc="Make figure #Moy", total=len(samples_list), ncols=75)) as progress_bar:
            Parallel(n_jobs= threads_number, verbose = 0)(delayed(Figures_CoCar.draw_figure)(sample, genes_list, out_directory + "/myCount", probes_dir) for sample in samples_list)
    
    if model:
        print("Step 5: myStatAnalysis...")
        if not os.path.exists(out_directory + "/myStatAnalysis"):
                os.mkdir(out_directory + "/myStatAnalysis")
        for gene in tqdm(genes_list, desc="Stat Analysis", total=len(genes_list), ncols=75):
            if not os.path.exists(out_directory + "/myStatAnalysis" + "/" + gene + "/Count"):
                os.makedirs(out_directory + "/myStatAnalysis" + "/" + gene + "/Count")
            path_list = glob.glob(out_directory + "/myCount" + "/[!#Moy!Undetermined]*/" + gene +"/Counts/" + "*.csv")
            for path in path_list:
                shutil.copyfile(path, out_directory + "/myStatAnalysis" + "/" + gene + "/Count" + "/" + path.split("/")[-1].split(".")[0] + ".count")
            subprocess.call("Rscript /mnt/recherche/RT-MLPA/RT-MLPA_scripts/CrossValAnalyzer.r -I " + out_directory + "/myStatAnalysis" + "/" + gene + "/Count/" + " -O " + out_directory + "/myStatAnalysis" + "/" + gene + "/" + gene + "_result.txt", stdout=subprocess.DEVNULL,  stderr=subprocess.STDOUT, shell = True)
    
    print("\n--------FINISH--------")
    print(f"--- Program {os.path.basename(__file__)} executed in {round((float(time.time() - start_time) / 60),2)} minutes ---")



if __name__ == "__main__":
    main()