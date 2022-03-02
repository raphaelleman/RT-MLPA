import os
import csv
import math
import time
import shutil
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
    parser.add_argument('-T', '--threads', dest='threads', type=str, default=os.cpu_count(),
                        help='number of threads to use for paralleling (default maximum cpu threads)')
    return parser.parse_args()    


def Count(list):
    count = {}.fromkeys(set(list),0)
    for i in list:
        count[i] += 1
    return count


def count_analysis(sample, directory, reads_dir, probes_dir):
    CoCar_Genes = ['BRCA1','BRCA2','PALB2','RAD51C','RAD51D', 'TP53', 'PTEN', 'CDH1', 'MSH2', 'MSH6', 'MLH1', 'PMS2', 'EPCAM']
    repS = directory + os.path.sep + sample
    if not os.path.exists(repS):
        os.mkdir(repS)

    Bilan = []
    print(sample)
    Reads_Ok_Precent = 0
    # Lecture du fichier reads
    os.chdir(reads_dir)
    namfile = str(sample) + '_Reads.csv'
    SeqT = pd.read_csv(namfile, sep=';' , encoding='latin-1', header = None)

    #SeqT = SeqT[0:500000]

    SeqT.columns = ['Seq']
    nb = len(SeqT)
    Bilan.append(namfile)
    Bilan.append('Reads_Totaux_FastQ')
    Bilan.append(nb)

    # Extraction UMI
    SeqT['UMI'] = SeqT['Seq'].str[2:9]

    # recherche des sondes et barcodes
    SeqT['Left'] = '_'
    SeqT['Right'] = '_'

    for Analysis in CoCar_Genes:
        repG = repS + os.path.sep + Analysis
        if not os.path.exists(repG):
            os.mkdir(repG)
            os.mkdir(repG + os.path.sep + 'Aligns')
            os.mkdir(repG + os.path.sep + 'Counts')
            os.mkdir(repG + os.path.sep + 'Figures')
            os.mkdir(repG + os.path.sep + 'Figures/Matrices')
        Bilan.append(Analysis)
        print(str(len(SeqT)) + ' Total Reads')
        print(Analysis)
        SeqT['Seq2'] = '_'
        SeqT['tmp'] = SeqT['Seq'].str[9:19]

        ## Chargement et tri des sondes (D/G)
        os.chdir(probes_dir)
        probes = Analysis + '_Probes.csv'
        sondes = pd.read_csv(probes, sep=';' , encoding='latin-1')
        sondes.index = sondes['Sonde']
        sondes['DG'] = '_'
        sondes['10bp'] = '_'
        sondes['S'] = '0'
        for i in range(0,len(sondes)):
            sondes['DG'][i] = (sondes['Sonde'][i])[len(sondes['Sonde'][i])-1]
            sondes['10bp'][i] = (sondes['Seq'][i])[0:10]
            sondes['S'][i] = len(sondes['Seq'][i])
        sondesG = sondes[sondes['DG'] == 'G']
        sondesD = sondes[sondes['DG'] == 'D']

        ## Alignement des sondes _ sauvegarde des fichiers d'alignement
        # Alignement _____________________________________________
        for sd in sondesG['Sonde']:
            sd = sd.replace(' ','')
            seq = str(sondesG['10bp'][sd])
            SeqT['Left'][SeqT['tmp'] == seq] = sd
            L1 = sondesG['S'][sd]
            SeqT['Seq2'][SeqT['Left'] == sd] = SeqT['Seq'].str[L1:]
        print('Left Probes ok')
        Bilan.append('Left_ok')
        Bilan.append(len(SeqT[SeqT['Left'] != '_']))

        SeqT['tmp'] = SeqT['Seq2'].str[9:19]
        for sd in sondesD['Sonde']:
            sd = sd.replace(' ','')
            seq = str(sondesD['10bp'][sd])
            SeqT['Right'][SeqT['tmp'] == seq] = sd
            L2 = sondesD['S'][sd] + 28
        print('Right Probes ok')
        Bilan.append('Right_ok')
        Bilan.append(len(SeqT[SeqT['Right'] != '_']))

        del SeqT['tmp']
        del SeqT['Seq2']

        # Sauvegarde
        os.chdir(repG + os.path.sep + 'Aligns')
        namfile = str(sample) + '_Align_UnFiltred_' + Analysis + '.csv'
        #SeqT.to_csv(namfile, sep = ';', index = False)

        # Filtres : barcodes non ok et pb sondes
        SeqT_Filtred = SeqT[SeqT['Left'] != '_']

        SeqT_Left = SeqT_Filtred[SeqT_Filtred['Right'] == '_']

        SeqT_Filtred = SeqT_Filtred[SeqT_Filtred['Right'] != '_']
        SeqT_Filtred = SeqT_Filtred.reset_index(drop = True)

        Reads_Ok_Precent = Reads_Ok_Precent + len(SeqT_Filtred)/nb*100

        print(str(len(SeqT_Filtred)) + ' Reads Ok   ' + str(round(Reads_Ok_Precent,1)) + ' %')

        #retirer les reads identifiés
        SeqT = SeqT[SeqT['Left'] == '_']

        # Sauvegarde
        namfile = str(sample) + '_Align_Left_Only' + Analysis + '.csv'
        #SeqT_Left.to_csv(namfile, sep = ';', index = False)

        namfile = str(sample) + '_Align_' + Analysis + '.csv'
        SeqT_Filtred.to_csv(namfile, sep = ';', index = False)

    # Comptage _______________________________________________
        # Chargement du fichier
        os.chdir(repG + os.path.sep + 'Aligns')
        namfile = str(sample) + '_Align_' + Analysis + '.csv'
        #SeqT_Filtred = pd.read_csv(namfile, sep=';' , encoding='latin-1')
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
                # Comptage UMI
            FusUMI = SeqT_Filtred['FusUMI'].tolist()
            UMI = pd.DataFrame(list(Count(FusUMI).items()))
            UMI.columns = ['Fusion', 'Count_UMI']

            UMI['Fusion'] = UMI['Fusion'].str[:-8]

            UMIcpt = UMI['Fusion'].tolist()
            cptUMI = pd.DataFrame(list(Count(UMIcpt).items()))
            cptUMI.columns = ['Fusion', 'Count_UMI']

                # Comptage UMI corrige : premier quartile 0.25 (TBOne = 0.20)
            quant = np.quantile(UMI['Count_UMI'],0.25)
            UMIcorr = UMI[UMI['Count_UMI'] >= quant]
            UMIcorr=UMIcorr.reset_index(drop = True)
            UMIcorrcpt = UMIcorr['Fusion'].tolist()
            cptUMIcorr = pd.DataFrame(list(Count(UMIcorrcpt).items()))
            cptUMIcorr.columns = ['Fusion', 'Count_UMI_Adj']

                # Construction tableau de comptage
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

            # Agglomerer les sondes degenerees
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

            namfile = str(sample) + '_Counts_' + Analysis + '.csv'
            os.chdir(repG + os.path.sep + 'Counts')
            Probes.to_csv(namfile, sep = ';')
            print('stop')

    os.chdir(repS)
    namfile = str(sample) + '_UnAligned.csv'
    #SeqT.to_csv(namfile, sep = ';', index = False)
    bilan = pd.DataFrame(Bilan)
    namfile = str(sample) + '_info.csv'
    bilan.to_csv(namfile, sep = ';', index = False)



def main():
    start_time = time.time()

    #get arguments
    args = get_arguments()
    directory = os.path.realpath(args.output) + os.path.sep + args.runName
    probes_dir = os.path.realpath(args.probes)
    reads_dir = os.path.realpath(args.input)
    threads_number = args.threads
    
    #create directory if it not exists
    if not os.path.exists(directory):
        os.mkdir(directory)

    #get sample list
    reads_list = [dir for dir in os.listdir(reads_dir)]
    sample_list = [read.split("_")[0] for read in reads_list]

    #launch analysis
    Parallel(n_jobs= threads_number, verbose = 1)(delayed(count_analysis)(sample, directory, reads_dir, probes_dir) for sample in sample_list)

    print("--- Program executed in {:2f} minutes ---".format(float(time.time() - start_time) / 60))



if __name__ == "__main__":
    main()