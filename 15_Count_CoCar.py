import os
import time
import glob
import argparse
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
pd.set_option('mode.chained_assignment', None)


def get_arguments():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('-R', '--runName', dest='runName', type=str,
                        help='Name of run')
    parser.add_argument('-O', '--output', dest='output', type=str,
                        help='/path/to/output')
    parser.add_argument('-P', '--probes', dest='probes', type=str,
                        help='/path/to/probes ref')
    parser.add_argument('-I', '--input', dest='input', type=str,
                        help='/path/to/XXXX_Reads.csv files')
    parser.add_argument('-G', '--genes_file', dest='gene_file', type=str,
                        help='/path/to/genes file')
    parser.add_argument('-T', '--threads', dest='threads', type=str, default=os.cpu_count(),
                        help='number of threads to use for paralleling (default maximum cpu threads)')
    return parser.parse_args()    


def Count(list):
    count = {}.fromkeys(set(list),0)
    for i in list:
        count[i] += 1
    return count


def count_analysis(sample, genes_list, directory, reads_dir, probes_dir):
    sample_dir = directory + os.path.sep + sample
    if not os.path.exists(sample_dir):
        os.mkdir(sample_dir)

    report = []
    print(sample)
    Reads_Ok_Precent = 0
    #Read file
    os.chdir(reads_dir)
    filename = str(sample) + '_Reads.csv'
    SeqT = pd.read_csv(filename, sep=';' , encoding='latin-1', header = None)

    #SeqT = SeqT[0:500000]

    SeqT.columns = ['Seq']
    nb = len(SeqT)
    report.append(filename)
    report.append('Reads_Totaux_FastQ')
    report.append(nb)

    #Extract UMI
    SeqT['UMI'] = SeqT['Seq'].str[2:9]

    #Search probes and barcodes
    SeqT['Left'] = '_'
    SeqT['Right'] = '_'

    for gene in genes_list:
        gene_dir = sample_dir + os.path.sep + gene
        if not os.path.exists(gene_dir):
            os.mkdir(gene_dir)
            os.mkdir(gene_dir + os.path.sep + 'Aligns')
            os.mkdir(gene_dir + os.path.sep + 'Counts')
            os.mkdir(gene_dir + os.path.sep + 'Figures')
            os.mkdir(gene_dir + os.path.sep + 'Figures/Matrices')
        report.append(gene)
        print(str(len(SeqT)) + ' Total Reads')
        print(gene)
        SeqT['Seq2'] = '_'
        SeqT['tmp'] = SeqT['Seq'].str[9:19]

        ##Load and sort probes (L/R)
        os.chdir(probes_dir)
        probe_name = gene + '_Probes.csv'
        probes = pd.read_csv(probe_name, sep=';' , encoding='latin-1')
        probes.index = probes['Sonde']
        probes['DG'] = '_'
        probes['10bp'] = '_'
        probes['S'] = '0'
        for i in range(0,len(probes)):
            probes['DG'][i] = (probes['Sonde'][i])[len(probes['Sonde'][i])-1]
            probes['10bp'][i] = (probes['Seq'][i])[0:10]
            probes['S'][i] = len(probes['Seq'][i])
        left_probes = probes[probes['DG'] == 'G']
        right_probes = probes[probes['DG'] == 'D']

        ## Align probes _ save alignment file
        # Align _____________________________________________
        for lp in left_probes['Sonde']:
            lp = lp.replace(' ','')
            seq = str(left_probes['10bp'][lp])
            SeqT['Left'][SeqT['tmp'] == seq] = lp
            L1 = left_probes['S'][lp]
            SeqT['Seq2'][SeqT['Left'] == lp] = SeqT['Seq'].str[L1:]
        print('Left Probes ok')
        report.append('Left_ok')
        report.append(len(SeqT[SeqT['Left'] != '_']))

        SeqT['tmp'] = SeqT['Seq2'].str[9:19]
        for rp in right_probes['Sonde']:
            rp = rp.replace(' ','')
            seq = str(right_probes['10bp'][rp])
            SeqT['Right'][SeqT['tmp'] == seq] = rp
            L2 = right_probes['S'][rp] + 28
        print('Right Probes ok')
        report.append('Right_ok')
        report.append(len(SeqT[SeqT['Right'] != '_']))

        del SeqT['tmp']
        del SeqT['Seq2']

        # Save
        os.chdir(gene_dir + os.path.sep + 'Aligns')
        filename = str(sample) + '_Align_UnFiltred_' + gene + '.csv'
        #SeqT.to_csv(filename, sep = ';', index = False)

        # Filtres : barcodes non ok et pb sondes
        SeqT_Filtred = SeqT[SeqT['Left'] != '_']

        SeqT_Left = SeqT_Filtred[SeqT_Filtred['Right'] == '_']

        SeqT_Filtred = SeqT_Filtred[SeqT_Filtred['Right'] != '_']
        SeqT_Filtred = SeqT_Filtred.reset_index(drop = True)

        Reads_Ok_Precent = Reads_Ok_Precent + len(SeqT_Filtred)/nb*100

        print(str(len(SeqT_Filtred)) + ' Reads Ok   ' + str(round(Reads_Ok_Precent,1)) + ' %')

        #Remove identified reads
        SeqT = SeqT[SeqT['Left'] == '_']

        #Save
        filename = str(sample) + '_Align_Left_Only' + gene + '.csv'
        #SeqT_Left.to_csv(filename, sep = ';', index = False)

        filename = str(sample) + '_Align_' + gene + '.csv'
        SeqT_Filtred.to_csv(filename, sep = ';', index = False)

    #Count _______________________________________________
        #Load file
        os.chdir(gene_dir + os.path.sep + 'Aligns')
        filename = str(sample) + '_Align_' + gene + '.csv'
        #SeqT_Filtred = pd.read_csv(filename, sep=';' , encoding='latin-1')
        del SeqT_Filtred['Seq']

        # Fusion des champs
        SeqT_Filtred['FusProbes'] = '_'
        SeqT_Filtred['FusUMI'] = '_'
        SeqT_Filtred['FusProbes'] = SeqT_Filtred['Left'] + '-' + SeqT_Filtred['Right']
        SeqT_Filtred['FusUMI'] = SeqT_Filtred['Left'] + '-' + SeqT_Filtred['Right'] + '-' + SeqT_Filtred['UMI']
        del SeqT_Filtred['UMI']
        del SeqT_Filtred['Left']
        del SeqT_Filtred['Right']

            # Comptage Full
        FusProbes = SeqT_Filtred['FusProbes'].tolist()
        Probes = pd.DataFrame(list(Count(FusProbes).items()))

        if Probes.empty == False :
            Probes.columns = ['Fusion', 'Count_Full']
                #Count UMI
            FusUMI = SeqT_Filtred['FusUMI'].tolist()
            UMI = pd.DataFrame(list(Count(FusUMI).items()))
            UMI.columns = ['Fusion', 'Count_UMI']

            UMI['Fusion'] = UMI['Fusion'].str[:-8]

            UMIcpt = UMI['Fusion'].tolist()
            cptUMI = pd.DataFrame(list(Count(UMIcpt).items()))
            cptUMI.columns = ['Fusion', 'Count_UMI']

                #Count corrected UMI : first quartile 0.25 (TBOne = 0.20)
            quant = np.quantile(UMI['Count_UMI'],0.25)
            UMIcorr = UMI[UMI['Count_UMI'] >= quant]
            UMIcorr=UMIcorr.reset_index(drop = True)
            UMIcorrcpt = UMIcorr['Fusion'].tolist()
            cptUMIcorr = pd.DataFrame(list(Count(UMIcorrcpt).items()))
            cptUMIcorr.columns = ['Fusion', 'Count_UMI_Adj']

                #Create count table
            Probes=Probes.reset_index(drop = True)
            Probes.index = Probes['Fusion']
            fusion = Probes['Fusion'].tolist()
            del Probes['Fusion']

            Probes['Count_UMI'] = '0'
            for fus in fusion:
                if (fus in cptUMI['Fusion'].tolist()) == True:
                    val = cptUMI['Count_UMI'][cptUMI['Fusion'] == fus]
                    val = (val.tolist())[0]
                    Probes['Count_UMI'][fus] = val

            Probes['Adj_UMI'] = '0'

            for fus in fusion:
                if (fus in cptUMIcorr['Fusion'].tolist()) == True:
                    val = cptUMIcorr['Count_UMI_Adj'][cptUMIcorr['Fusion'] == fus]
                    val = (val.tolist())[0]
                    Probes['Adj_UMI'][fus] = val

            #Agglomerate degenerates probes
            for probes_2 in Probes.index:
                if '_2' in probes_2:
                    probes_3 = probes_2.replace('_2','')
                    if probes_2 in Probes.index and probes_3 in Probes.index:
                        Probes['Count_Full'][probes_3] = int(Probes['Count_Full'][probes_2]) + int(Probes['Count_Full'][probes_3])
                        Probes['Count_UMI'][probes_3] = int(Probes['Count_UMI'][probes_2]) + int(Probes['Count_UMI'][probes_3])
                        Probes['Adj_UMI'][probes_3] = int(Probes['Adj_UMI'][probes_2]) + int(Probes['Adj_UMI'][probes_3])
                    if probes_3 not in Probes.index:
                        Probes = Probes.rename(index={probes_2: probes_3})
                    if probes_2 in Probes.index:
                        Probes = Probes.drop(index = probes_2)

            filename = str(sample) + '_Counts_' + gene + '.csv'
            os.chdir(gene_dir + os.path.sep + 'Counts')
            Probes.to_csv(filename, sep = ';')
            print('stop')

    os.chdir(sample_dir)
    filename = str(sample) + '_UnAligned.csv'
    #SeqT.to_csv(filename, sep = ';', index = False)
    report = pd.DataFrame(report)
    filename = str(sample) + '_info.csv'
    report.to_csv(filename, sep = ';', index = False)



def main():
    start_time = time.time()

    #get arguments
    args = get_arguments()
    directory = os.path.realpath(args.output) + os.path.sep + args.runName
    probes_dir = os.path.realpath(args.probes)
    reads_dir = os.path.realpath(args.input)
    threads_number = args.threads
    gene_file = os.path.realpath(args.gene_file)
    
    #create directory if it not exists
    if not os.path.exists(directory):
        os.mkdir(directory)

    #get sample list
    os.chdir(reads_dir)
    reads_list = glob.glob("*.csv")
    samples_list = [read.split("_")[0] for read in reads_list]
    
    # get gene list from file
    genes_list = []
    with open(gene_file, "r") as file:
        for line in file:
            genes_list.append(line.strip())
    print(genes_list)
    

    #launch analysis
    Parallel(n_jobs= threads_number, verbose = 1)(delayed(count_analysis)(sample, genes_list, directory, reads_dir, probes_dir) for sample in samples_list)

    print(f"--- Program executed in {float(time.time() - start_time) / 60} minutes ---")



if __name__ == "__main__":
    main()