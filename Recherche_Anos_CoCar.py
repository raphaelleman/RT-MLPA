import os
import time
import argparse
import numpy as np
import pandas as pd
from statistics import mean
from joblib import Parallel, delayed
pd.set_option('mode.chained_assignment', None)
import warnings
warnings.filterwarnings('ignore')

def get_arguments():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('-R', '--runName', dest='runName', type=str,
                        help='Name of run')
    parser.add_argument('-O', '--output', dest='output', type=str,
                        help='/path/to/output')
    parser.add_argument('-P', '--probes', dest='probes', type=str,
                        help='/path/to/probes ref')
    parser.add_argument('-I', '--input', dest='input', type=str,
                        help='/path/to/count directory')
    parser.add_argument('-G', '--genes_file', dest='gene_file', type=str,
                        help='/path/to/genes file')
    parser.add_argument('-T', '--threads', dest='threads', type=str, default=os.cpu_count(),
                        help='number of threads to use for paralleling (default maximum cpu threads)')
    return parser.parse_args()  


def search_anos(gene, index_matrix, directory, run_name, counts_dir, probes_dir):
    #print('')
    #print(gene)
    #print('Creation de la matrice')

    os.chdir(probes_dir)
    col_matrix = []
    Probes = pd.read_csv(gene + '_Probes.csv', sep=';' , encoding='latin-1')
    Probes = Probes.Sonde
    left = []
    right = []
    for j in Probes:
        if j[-1] == 'G':
            left.append(j)
        if j[-1] == 'D':
            right.append(j)
    for j in left:
        for k in right:
            col_matrix.append(j + '_' + k)

    matrix = pd.DataFrame(columns = col_matrix, index = index_matrix)
    #print(matrix)
    # parcours des patients
    for j in index_matrix:
        Mat_dir = counts_dir + os.path.sep + str(j) + os.path.sep + str(gene) + os.path.sep + 'Figures/Matrices'

        if len(os.listdir(Mat_dir)) > 0:
            Mat_name = Mat_dir + os.path.sep + str(j) + '_Counts_' + str(gene) + '.csv'
            Mat = pd.read_csv(Mat_name, sep = ';', index_col = 'Sonde')

            # parcours des jonctions
            for k in Mat.index:
                for l in Mat.columns:
                    col = k + '_' + l
                    matrix[col][j] = Mat[l][k]

    matrix.index = matrix.index + '_' + run_name
    os.chdir(directory)

    matrix.insert(0,'Tot_Reads', 0)
    for k in matrix.index:
        coeff_norm = 0
        for l in matrix.columns[1:]:
            if matrix[l][k] == matrix[l][k]:
                coeff_norm = coeff_norm + matrix[l][k]
        matrix['Tot_Reads'][k] = coeff_norm

    matrix.to_csv('Matrice_Full_' + str(gene) + '_' + run_name + '.csv', sep = ';', index = True)

    ###### normalisation
    # Option 1 : par le total des reads du marqueur
    # matrix_2 = Matrice.copy()
    # for k in matrix_2.index:
    #     for l in matrix_2.columns[0:]:
    #         matrix_2[l][k] = matrix_2[l][k] / matrix_2['Tot_Reads'][k] * 1000000
    # matrix_2.to_csv('Matrice_Full_' + str(i) + '_Norm_' + run + '.csv', sep = ';', index = True)

    # Option 2 : par le nombre de reads des jonctions normales
    # creation de la liste des jonctions normales
    os.chdir(probes_dir)
    gene_file = pd.read_csv(gene + '_Design.csv', sep=';' , encoding='latin-1')
    listG = []
    listD = []
    for k in gene_file.columns:
        if k[0:len(gene)] == gene and k[len(gene)] == 'E' and k[len(gene)+1].isnumeric():
            if k[len(k)-1] == 'G' and k[len(k)-2].isnumeric():
                if '_' not in k and k != 'EPCAME9G':
                    listG.append(k)
            if k[len(k)-1] == 'D' and k[len(k)-2].isnumeric():
                if '_' not in k and k != 'GREM1E3D':
                    listD.append(k)
    del listG[len(listG)-1]
    del listD[0]
    ListJunctions = []
    for k in range(len(listG)):
        junc = listG[k] + '_' + listD[k]
        if junc != 'TP53E11G_TP53E12D' and junc != 'TP53E12G_TP53E13D' and junc != 'TP53E13G_TP53E14D' and junc != 'TP53E14G_TP53E15D' and listG[k] != 'MSH2E16G' and listG[k] != 'TP53E11G':
            ListJunctions.append(junc)

    # normalisation
    #print('Normalisation')
    matrix_2 = matrix.copy()
    for k in matrix_2.index:
        n = 0
        for l in ListJunctions:
            n = n + matrix_2[l][k]
        for l in matrix_2.columns[0:]:
            matrix_2[l][k] = matrix_2[l][k] / n * 1000000

    os.chdir(directory)
    matrix_2.to_csv('Matrice_Full_' + str(gene) + '_Norm_' + run_name + '.csv', sep = ';', index = True)

    ###### Creation du profil de reference
    #print('Generation du profil de reference')
    for Case in matrix_2.index:
        patient =  Case
        #print(patient)
        ExJ = []
        Counts = []
        for Junction in matrix.columns[1:]:
            #print(Junction)
            #Moy = mean(Matrice[Junction])
            Moy = matrix_2[Junction].quantile([0.50]).values[0]
            #print(Moy)
            if Moy != 0:
                ExJ.append(Junction)
                Counts.append(Moy)
    Matrice_Moy = pd.DataFrame()
    Matrice_Moy['Fusion'] = ExJ
    Matrice_Moy['Count_Full'] = Counts

    os.chdir(directory + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Counts')
    Matrice_Moy.to_csv('#Moy_Counts_' + gene + '.csv', sep = ';', index = False, header = True)

    ######## recherche des anomalies de splice
    #print('Detection des anomalies de splice')
    AnoSpliceList = []

    # 1 : Variations brutales
    # parcours des Cas :
    for Case in matrix_2.index:

        # Niveau moyen d'expression des joonctions normales
        ExpMoyenne = 0
        for junc in ListJunctions:
            ExpMoyenne = ExpMoyenne + matrix_2[junc][Case]
        ExpMoyenne = ExpMoyenne / len(ListJunctions)

        # 1a : test outliers isolés low ou high
            # low q25 - 1.5 IQR et < 80%
            # high q75 + 1.5 IQR et > 150%
        nbAno = False
        # parcours des jonctions :
        l = 0
        for k in ListJunctions:
            SimpleK = k.replace(gene,'')
            SimpleK = SimpleK.replace('G','')
            SimpleK = SimpleK.replace('D','')
            SimpleK = SimpleK.replace('_','-')

            Val = matrix_2[ListJunctions[l]][Case]
            q25 = matrix_2[ListJunctions[l]].quantile([0.25]).values[0]
            q50 = matrix_2[ListJunctions[l]].quantile([0.50]).values[0]
            q75 = matrix_2[ListJunctions[l]].quantile([0.75]).values[0]
            l += 1

            # low
            if Val < q25-(1.5*(q75 - q25)) and Val/q50*100 < 80 :
                Line = 'Low Outlier ' + SimpleK + ' [ ' + str(int(Val)) + ' vs ( ' + str(int(q25)) + '-' + str(int(q50)) + '-' + str(int(q75)) + ' ) / ' + str(int(Val/q50*100)) + ' % ]'
                #Line = 'Low Outlier ' + SimpleK + ' [ ' + str(int(Val/q50*100)) + ' % ]'

                if Val/q50*100 < 50:
                    Line = Line + ' *****'

                #print(Case[0:4])
                #print(Line)
                if nbAno == False :
                    AnoSpliceList.append(Case)
                    nbAno = True
                AnoSpliceList.append(Line)

            # high
            if Val > q75+(1.5*(q75 - q25)) and Val/q50*100 > 150 and q50 >0 :
                Line = 'High Outlier ' + SimpleK + ' [ ' + str(int(Val)) + ' vs ( ' + str(int(q25)) + '-' + str(int(q50)) + '-' + str(int(q75)) + ' ) / ' + str(int(Val/q50*100)) + ' % ]'
                #Line = 'High Outlier ' + SimpleK + ' [ ' + str(int(Val/q50*100)) + ' % ]'

                if Val/q50*100 > 300:
                    Line = Line + ' *****'

                #print(Case[0:4])
                #print(Line)
                if nbAno == False :
                    AnoSpliceList.append(Case)
                    nbAno = True
                AnoSpliceList.append(Line)

        # 1b : test exon skipping
        # parcours des jonctions :
        l = 0
        for SD_gauche in listG:
            for SD_droite in listD:
                Tested_Junc = SD_gauche + '_' + SD_droite
                if Tested_Junc != 'TP53E11G_TP53E12D' and Tested_Junc != 'TP53E12G_TP53E13D' and Tested_Junc != 'TP53E13G_TP53E14D' and Tested_Junc != 'TP53E14G_TP53E15D' and SD_gauche != 'MSH2E16G'  and SD_gauche != 'TP53E11G' and SD_gauche != 'TP53E13G' and SD_gauche != 'TP53E14G':

                    SimpleK_2 = Tested_Junc
                    SimpleK_2 = SimpleK_2.replace(gene,'')
                    SimpleK_2 = SimpleK_2.replace('G','')
                    SimpleK_2 = SimpleK_2.replace('D','')
                    SimpleK_2 = SimpleK_2.replace('_','-')

                    if listD.index(SD_droite) > listG.index(SD_gauche):
                        event = 'ExSkipping '
                    else:
                        event = 'Backsplice '

                    Val_2 = matrix_2[Tested_Junc][Case]
                    q25_2 = matrix_2[Tested_Junc].quantile([0.25]).values[0]
                    q50_2 = matrix_2[Tested_Junc].quantile([0.50]).values[0]
                    q75_2 = matrix_2[Tested_Junc].quantile([0.75]).values[0]

                    Reads_2 = matrix_2['Tot_Reads'][Case]

                    if Val_2 > (q75_2 + (1.5*(q75_2 - q25_2))) and event == 'ExSkipping ' and (Val_2/q50_2*100 > 500 or q50_2 == 0):
                        # skipping : au moins 5% du niveau d'expression moyen des jonctions normales
                        if Val_2 > Reads_2 / (len(ListJunctions) * 20):
                            if q50_2 >0 :
                                Line = event + SimpleK_2 + ' [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / ' + str(int(Val_2/q50_2*100)) + ' % ]'
                            else:
                                Line = event + SimpleK_2 + ' [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / inf % ]'

                            if event == 'ExSkipping ':
                                Line = Line + ' *****'

                            #print(Case[0:4])
                            #print(Line)
                            if nbAno == False :
                                AnoSpliceList.append(Case)
                                nbAno = True
                            AnoSpliceList.append(Line)

        # 2a : test outliers moderes <q25 mais contigus et < 70%
        # parcours des jonctions normales contigues:
        l = 0
        for k in ListJunctions[0:len(ListJunctions)-1]:
            SimpleK = k.replace(gene,'')
            SimpleK = SimpleK.replace('G','')
            SimpleK = SimpleK.replace('D','')
            SimpleK = SimpleK.replace('_','-')

            SimpleK_2 = ListJunctions[l+1]
            SimpleK_2 = SimpleK_2.replace(gene,'')
            SimpleK_2 = SimpleK_2.replace('G','')
            SimpleK_2 = SimpleK_2.replace('D','')
            SimpleK_2 = SimpleK_2.replace('_','-')

            Val_1 = matrix_2[ListJunctions[l]][Case]
            q25_1 = matrix_2[ListJunctions[l]].quantile([0.25]).values[0]
            q50_1 = matrix_2[ListJunctions[l]].quantile([0.50]).values[0]
            q75_1 = matrix_2[ListJunctions[l]].quantile([0.75]).values[0]

            Val_2 = matrix_2[ListJunctions[l+1]][Case]
            q25_2 = matrix_2[ListJunctions[l+1]].quantile([0.25]).values[0]
            q50_2 = matrix_2[ListJunctions[l+1]].quantile([0.50]).values[0]
            q75_2 = matrix_2[ListJunctions[l+1]].quantile([0.75]).values[0]

            l += 1

            if Val_1 < q25_1 and Val_2 < q25_2 and Val_1/q50_1*100 < 70 and Val_2/q50_2*100 < 70 :

                Line = 'Low Contig ' + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1)) + ' vs ( ' + str(int(q25_1)) + '-' + str(int(q50_1)) + '-' + str(int(q75_1)) + ' ) / ' + str(int(Val_1/q50_1*100)) + ' % ] and [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / ' + str(int(Val_2/q50_2*100)) + ' % ]'
                #Line = 'Low Contig ' + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1/q50_1*100)) + ' % and ' + str(int(Val_2/q50_2*100)) + ' % ]'

                if Val_1/q50_1*100 < 50 or Val_2/q50_2*100 <50:
                    Line = Line + ' *****'

                #print(Case[0:4])
                #print(Line)
                if nbAno == False :
                    AnoSpliceList.append(Case)
                    nbAno = True
                AnoSpliceList.append(Line)

        # 2b : test outliers moderes <q25 mais associees à des abberations > q75 + 1.5 IQR et significative (val > nbr_reads jonctions normales / nbre jonctions normales / 100)
        # parcours des jonctions :
        l = 0
        for k in ListJunctions:
            SimpleK = k.replace(gene,'')
            SimpleK = SimpleK.replace('G','')
            SimpleK = SimpleK.replace('D','')
            SimpleK = SimpleK.replace('_','-')

            # sonde de gauche
            SD_gauche = k[0:k.find('_')]

            Val_1 = matrix_2[ListJunctions[l]][Case]
            q25_1 = matrix_2[ListJunctions[l]].quantile([0.25]).values[0]
            q50_1 = matrix_2[ListJunctions[l]].quantile([0.50]).values[0]
            q75_1 = matrix_2[ListJunctions[l]].quantile([0.75]).values[0]

            for SD_droite in listD:
                Tested_Junc = SD_gauche + '_' + SD_droite

                SimpleK_2 = Tested_Junc
                SimpleK_2 = SimpleK_2.replace(gene,'')
                SimpleK_2 = SimpleK_2.replace('G','')
                SimpleK_2 = SimpleK_2.replace('D','')
                SimpleK_2 = SimpleK_2.replace('_','-')

                if listD.index(SD_droite) > listG.index(SD_gauche):
                    event = 'ExSkipping '
                else:
                    event = 'Backsplice '

                Val_2 = matrix_2[Tested_Junc][Case]
                q25_2 = matrix_2[Tested_Junc].quantile([0.25]).values[0]
                q50_2 = matrix_2[Tested_Junc].quantile([0.50]).values[0]
                q75_2 = matrix_2[Tested_Junc].quantile([0.75]).values[0]

                Reads_2 = matrix_2['Tot_Reads'][Case]

                if Val_1 < q25_1 and Val_2 > (q75_2 + (1.5*(q75_2 - q25_2))):
                    # Selection skipping ou backsplice : au moins 10% du niveau d'expression moyen des jonctions normales
                    if Val_2 > Reads_2 / (len(ListJunctions) * 10):
                        if q50_2 >0 :
                            Line = 'Low + ' + event + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1)) + ' vs ( ' + str(int(q25_1)) + '-' + str(int(q50_1)) + '-' + str(int(q75_1)) + ' ) / ' + str(int(Val_1/q50_1*100)) + ' % ] and [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / ' + str(int(Val_2/q50_2*100)) + ' % ]'
                            #Line = 'Low + ' + event + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1/q50_1*100)) + ' % and ' + str(int(Val_2/q50_2*100)) + ' % ]'
                        else:
                            Line = 'Low + ' + event + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1)) + ' vs ( ' + str(int(q25_1)) + '-' + str(int(q50_1)) + '-' + str(int(q75_1)) + ' ) / ' + str(int(Val_1/q50_1*100)) + ' % ] and [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / inf % ]'
                            #Line = 'Low + ' + event + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1/q50_1*100)) + ' % and inf % ]'

                        if event == 'ExSkipping ':
                            Line = Line + ' *****'

                        #print(Case[0:4])
                        #print(Line)
                        if nbAno == False :
                            AnoSpliceList.append(Case)
                            nbAno = True
                        AnoSpliceList.append(Line)

            l += 1

        if nbAno == True:
            AnoSpliceList.append('')

    os.chdir(directory)
    AnoSpliceList = pd.DataFrame(AnoSpliceList)
    AnoSpliceList.to_csv('Ano_' + gene + '.csv', sep = ';', header = False, index = False)



def main():
    start_time = time.time()

    #get arguments
    args = get_arguments()
    run_name = args.runName
    directory = os.path.realpath(args.output) + os.path.sep + run_name
    probes_dir = os.path.realpath(args.probes)
    counts_dir = os.path.realpath(args.input)
    gene_file = os.path.realpath(args.gene_file)
    threads_number = args.threads

     # get gene list from file
    genes_list = []
    with open(gene_file, "r") as file:
        for line in file:
            genes_list.append(line.strip())
    print(genes_list)

    #create directories
    if not os.path.exists(directory):
        os.mkdir(directory)
    os.mkdir(directory + os.path.sep + '#Moy')
    for gene in genes_list:
        os.mkdir(directory + os.path.sep + '#Moy' + os.path.sep + gene)
        os.mkdir(directory + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Counts')
        os.mkdir(directory + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Figures')
        os.mkdir(directory + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices')

    #get sample list
    index_matrix = []
    for file in os.listdir(counts_dir):
        if file[0:7] != 'Matrice' and file[0:4] != 'Ano_' and file[0:4] != '#Moy':
            index_matrix.append(file)

    #launch analysis
    Parallel(n_jobs= threads_number, verbose = 1)(delayed(search_anos)(gene, index_matrix, directory, run_name, counts_dir, probes_dir) for gene in genes_list)

    print(f"--- Program executed in {float(time.time() - start_time) / 60} minutes ---")



if __name__ == "__main__":
    main()