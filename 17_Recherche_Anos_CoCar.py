import os
import shutil
import csv
import math
import numpy as np
import pandas as pd
import time
from statistics import mean
import argparse
pd.set_option('mode.chained_assignment', None)

#run = '2021_04_14_CoCar1'
#run = '2021_06_04_CoCar2'
#run = '2021_09_15_CoCar3'
my_parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

my_parser.add_argument('-R','--runName',type=str,
                       help='Name of run')

my_parser.add_argument('-O','--output',type=str,
                       help='/path/to/output')

my_parser.add_argument('-P','--probes',type=str,
                       help='/path/to/probes ref')

my_parser.add_argument('-I','--input',type=str,
                       help='/path/to/count directory')

my_parser.add_argument('-v',
                       '--verbose',
                       action='store_true',
                       help='an optional argument')

args = my_parser.parse_args()
#run = '2021_04_14_CoCar1'
#run = '2021_06_04_CoCar2'
#run = '2021_09_15_CoCar3'
#run = '2021_09_27_CoCar4'
run = args.runName

repCount = os.path.realpath(args.input)

#CoCar_Genes = ['APC','AXIN2','BMPR1A','BUB1','EPCAM','FAN1','GALNT12','GREM1','MLH1','MSH2','MSH3','MSH6','MUTYH','NTHL1','PMS2','POLD1','POLE','PTEN','RNF43','RPS20','SMAD4','STK11','TP53']
#CoCar_Genes = ['APC','AXIN2','BMPR1A','BUB1','EPCAM','FAN1','GALNT12','GREM1','MLH1','MSH2','MSH3','MSH6','MUTYH','PMS2','POLD1','PTEN','RNF43','RPS20','SMAD4','STK11','TP53']
#CoCar_Genes = ['APC','BMPR1A','EPCAM','MLH1','MSH2','MSH6','MUTYH','NTHL1','PMS2','POLD1','POLE','PTEN','SMAD4','STK11','TP53','CDH1','PALB2','RAD51C','RAD51D','BRCA1','BRCA2']

CoCar_Genes = ['BRCA1','BRCA2','PALB2','RAD51C','RAD51D']

### Creation de la matrice globale
rep = os.path.realpath(args.output)
rep = rep+ os.path.sep + run
repProbes = os.path.realpath(args.probes)

### Creation des repertoire de la figure de reference
if not os.path.exists(rep):
    os.mkdir(rep)
    os.mkdir(rep + os.path.sep + '#Moy')
    for gene in CoCar_Genes:
        os.mkdir(rep + os.path.sep + '#Moy' + os.path.sep + gene)
        os.mkdir(rep + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Counts')
        os.mkdir(rep + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Figures')
        os.mkdir(rep + os.path.sep + '#Moy' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices')

# indexes : echantillons
Matrice_indexes = []
for file in os.listdir(repCount):
    if file[0:7] != 'Matrice' and file[0:4] != 'Ano_' and file[0:4] != '#Moy':
        Matrice_indexes.append(file)

# colonnes : ligations
for i in CoCar_Genes:
    print('')
    print(i)
    print('Creation de la matrice')

    os.chdir(repProbes)
    Matrice_colonnes = []
    Probes = pd.read_csv(i + '_Probes.csv', sep=';' , encoding='latin-1')
    Probes = Probes.Sonde
    Gauche = []
    Droite = []
    for j in Probes:
        if j[-1] == 'G':
            Gauche.append(j)
        if j[-1] == 'D':
            Droite.append(j)
    for j in Gauche:
        for k in Droite:
            Matrice_colonnes.append(j + '_' + k)

    Matrice = pd.DataFrame(columns = Matrice_colonnes, index = Matrice_indexes)
    # parcours des patients
    for j in Matrice_indexes:
        Mat_dir = repCount + os.path.sep + str(j) + os.path.sep + str(i) + os.path.sep + 'Figures/Matrices'

        if len(os.listdir(Mat_dir)) > 0:
            Mat_name = Mat_dir + os.path.sep + str(j) + '_Counts_' + str(i) + '.csv'
            Mat = pd.read_csv(Mat_name, sep = ';', index_col = 'Sonde')

            # parcours des jonctions
            for k in Mat.index:
                for l in Mat.columns:
                    col = k + '_' + l
                    Matrice[col][j] = Mat[l][k]

    Matrice.index = Matrice.index + '_' + run
    os.chdir(rep)

    Matrice.insert(0,'Tot_Reads', 0)
    for k in Matrice.index:
        coeff_norm = 0
        for l in Matrice.columns[1:]:
            if Matrice[l][k] == Matrice[l][k]:
                coeff_norm = coeff_norm + Matrice[l][k]
        Matrice['Tot_Reads'][k] = coeff_norm

    Matrice.to_csv('Matrice_Full_' + str(i) + '_' + run + '.csv', sep = ';', index = True)

    ###### normalisation
    # Option 1 : par le total des reads du marqueur
    # Matrice2 = Matrice.copy()
    # for k in Matrice2.index:
    #     for l in Matrice2.columns[0:]:
    #         Matrice2[l][k] = Matrice2[l][k] / Matrice2['Tot_Reads'][k] * 1000000
    # Matrice2.to_csv('Matrice_Full_' + str(i) + '_Norm_' + run + '.csv', sep = ';', index = True)

    # Option 2 : par le nombre de reads des jonctions normales
    # creation de la liste des jonctions normales
    os.chdir(repProbes)
    gene = pd.read_csv(i + '_Design.csv', sep=';' , encoding='latin-1')
    listG = []
    listD = []
    for k in gene.columns:
        if k[0:len(i)] == i and k[len(i)] == 'E' and k[len(i)+1].isnumeric():
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
    print('Normalisation')
    Matrice2 = Matrice.copy()
    for k in Matrice2.index:
        n = 0
        for l in ListJunctions:
            n = n + Matrice2[l][k]
        for l in Matrice2.columns[0:]:
            Matrice2[l][k] = Matrice2[l][k] / n * 1000000

    os.chdir(rep)
    Matrice2.to_csv('Matrice_Full_' + str(i) + '_Norm_' + run + '.csv', sep = ';', index = True)

    ###### Creation du profil de reference
    print('Generation du profil de reference')
    for Case in Matrice2.index:
        patient =  Case[0:4]
        #print(patient)
        ExJ = []
        Counts = []
        for Junction in Matrice.columns[1:]:
            #print(Junction)
            #Moy = mean(Matrice[Junction])
            Moy = Matrice2[Junction].quantile([0.50]).values[0]
            #print(Moy)
            if Moy != 0:
                ExJ.append(Junction)
                Counts.append(Moy)
    Matrice_Moy = pd.DataFrame()
    Matrice_Moy['Fusion'] = ExJ
    Matrice_Moy['Count_Full'] = Counts

    os.chdir(rep + os.path.sep + '#Moy' + os.path.sep + i + os.path.sep + 'Counts')
    Matrice_Moy.to_csv('#Moy_Counts_' + i + '.csv', sep = ';', index = False, header = True)

    ######## recherche des anomalies de splice
    print('Detection des anomalies de splice')
    AnoSpliceList = []

    # 1 : Variations brutales
    # parcours des Cas :
    for Case in Matrice2.index:

        # Niveau moyen d'expression des joonctions normales
        ExpMoyenne = 0
        for junc in ListJunctions:
            ExpMoyenne = ExpMoyenne + Matrice2[junc][Case]
        ExpMoyenne = ExpMoyenne / len(ListJunctions)

        # 1a : test outliers isolés low ou high
            # low q25 - 1.5 IQR et < 80%
            # high q75 + 1.5 IQR et > 150%
        nbAno = False
        # parcours des jonctions :
        l = 0
        for k in ListJunctions:
            SimpleK = k.replace(i,'')
            SimpleK = SimpleK.replace('G','')
            SimpleK = SimpleK.replace('D','')
            SimpleK = SimpleK.replace('_','-')

            Val = Matrice2[ListJunctions[l]][Case]
            q25 = Matrice2[ListJunctions[l]].quantile([0.25]).values[0]
            q50 = Matrice2[ListJunctions[l]].quantile([0.50]).values[0]
            q75 = Matrice2[ListJunctions[l]].quantile([0.75]).values[0]
            l += 1

            # low
            if Val < q25-(1.5*(q75 - q25)) and Val/q50*100 < 80 :
                Line = 'Low Outlier ' + SimpleK + ' [ ' + str(int(Val)) + ' vs ( ' + str(int(q25)) + '-' + str(int(q50)) + '-' + str(int(q75)) + ' ) / ' + str(int(Val/q50*100)) + ' % ]'
                #Line = 'Low Outlier ' + SimpleK + ' [ ' + str(int(Val/q50*100)) + ' % ]'

                if Val/q50*100 < 50:
                    Line = Line + ' *****'

                print(Case[0:4])
                print(Line)
                if nbAno == False :
                    AnoSpliceList.append(Case[0:4])
                    nbAno = True
                AnoSpliceList.append(Line)

            # high
            if Val > q75+(1.5*(q75 - q25)) and Val/q50*100 > 150 and q50 >0 :
                Line = 'High Outlier ' + SimpleK + ' [ ' + str(int(Val)) + ' vs ( ' + str(int(q25)) + '-' + str(int(q50)) + '-' + str(int(q75)) + ' ) / ' + str(int(Val/q50*100)) + ' % ]'
                #Line = 'High Outlier ' + SimpleK + ' [ ' + str(int(Val/q50*100)) + ' % ]'

                if Val/q50*100 > 300:
                    Line = Line + ' *****'

                print(Case[0:4])
                print(Line)
                if nbAno == False :
                    AnoSpliceList.append(Case[0:4])
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
                    SimpleK_2 = SimpleK_2.replace(i,'')
                    SimpleK_2 = SimpleK_2.replace('G','')
                    SimpleK_2 = SimpleK_2.replace('D','')
                    SimpleK_2 = SimpleK_2.replace('_','-')

                    if listD.index(SD_droite) > listG.index(SD_gauche):
                        event = 'ExSkipping '
                    else:
                        event = 'Backsplice '

                    Val_2 = Matrice2[Tested_Junc][Case]
                    q25_2 = Matrice2[Tested_Junc].quantile([0.25]).values[0]
                    q50_2 = Matrice2[Tested_Junc].quantile([0.50]).values[0]
                    q75_2 = Matrice2[Tested_Junc].quantile([0.75]).values[0]

                    Reads_2 = Matrice2['Tot_Reads'][Case]

                    if Val_2 > (q75_2 + (1.5*(q75_2 - q25_2))) and event == 'ExSkipping ' and (Val_2/q50_2*100 > 500 or q50_2 == 0):
                        # skipping : au moins 5% du niveau d'expression moyen des jonctions normales
                        if Val_2 > Reads_2 / (len(ListJunctions) * 20):
                            if q50_2 >0 :
                                Line = event + SimpleK_2 + ' [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / ' + str(int(Val_2/q50_2*100)) + ' % ]'
                            else:
                                Line = event + SimpleK_2 + ' [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / inf % ]'

                            if event == 'ExSkipping ':
                                Line = Line + ' *****'

                            print(Case[0:4])
                            print(Line)
                            if nbAno == False :
                                AnoSpliceList.append(Case[0:4])
                                nbAno = True
                            AnoSpliceList.append(Line)

        # 2a : test outliers moderes <q25 mais contigus et < 70%
        # parcours des jonctions normales contigues:
        l = 0
        for k in ListJunctions[0:len(ListJunctions)-1]:
            SimpleK = k.replace(i,'')
            SimpleK = SimpleK.replace('G','')
            SimpleK = SimpleK.replace('D','')
            SimpleK = SimpleK.replace('_','-')

            SimpleK_2 = ListJunctions[l+1]
            SimpleK_2 = SimpleK_2.replace(i,'')
            SimpleK_2 = SimpleK_2.replace('G','')
            SimpleK_2 = SimpleK_2.replace('D','')
            SimpleK_2 = SimpleK_2.replace('_','-')

            Val_1 = Matrice2[ListJunctions[l]][Case]
            q25_1 = Matrice2[ListJunctions[l]].quantile([0.25]).values[0]
            q50_1 = Matrice2[ListJunctions[l]].quantile([0.50]).values[0]
            q75_1 = Matrice2[ListJunctions[l]].quantile([0.75]).values[0]

            Val_2 = Matrice2[ListJunctions[l+1]][Case]
            q25_2 = Matrice2[ListJunctions[l+1]].quantile([0.25]).values[0]
            q50_2 = Matrice2[ListJunctions[l+1]].quantile([0.50]).values[0]
            q75_2 = Matrice2[ListJunctions[l+1]].quantile([0.75]).values[0]

            l += 1

            if Val_1 < q25_1 and Val_2 < q25_2 and Val_1/q50_1*100 < 70 and Val_2/q50_2*100 < 70 :

                Line = 'Low Contig ' + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1)) + ' vs ( ' + str(int(q25_1)) + '-' + str(int(q50_1)) + '-' + str(int(q75_1)) + ' ) / ' + str(int(Val_1/q50_1*100)) + ' % ] and [ ' + str(int(Val_2)) + ' vs ( ' + str(int(q25_2)) + '-' + str(int(q50_2)) + '-' + str(int(q75_2)) + ' ) / ' + str(int(Val_2/q50_2*100)) + ' % ]'
                #Line = 'Low Contig ' + SimpleK + ' + ' + SimpleK_2 + ' [ ' + str(int(Val_1/q50_1*100)) + ' % and ' + str(int(Val_2/q50_2*100)) + ' % ]'

                if Val_1/q50_1*100 < 50 or Val_2/q50_2*100 <50:
                    Line = Line + ' *****'

                print(Case[0:4])
                print(Line)
                if nbAno == False :
                    AnoSpliceList.append(Case[0:4])
                    nbAno = True
                AnoSpliceList.append(Line)

        # 2b : test outliers moderes <q25 mais associees à des abberations > q75 + 1.5 IQR et significative (val > nbr_reads jonctions normales / nbre jonctions normales / 100)
        # parcours des jonctions :
        l = 0
        for k in ListJunctions:
            SimpleK = k.replace(i,'')
            SimpleK = SimpleK.replace('G','')
            SimpleK = SimpleK.replace('D','')
            SimpleK = SimpleK.replace('_','-')

            # sonde de gauche
            SD_gauche = k[0:k.find('_')]

            Val_1 = Matrice2[ListJunctions[l]][Case]
            q25_1 = Matrice2[ListJunctions[l]].quantile([0.25]).values[0]
            q50_1 = Matrice2[ListJunctions[l]].quantile([0.50]).values[0]
            q75_1 = Matrice2[ListJunctions[l]].quantile([0.75]).values[0]

            for SD_droite in listD:
                Tested_Junc = SD_gauche + '_' + SD_droite

                SimpleK_2 = Tested_Junc
                SimpleK_2 = SimpleK_2.replace(i,'')
                SimpleK_2 = SimpleK_2.replace('G','')
                SimpleK_2 = SimpleK_2.replace('D','')
                SimpleK_2 = SimpleK_2.replace('_','-')

                if listD.index(SD_droite) > listG.index(SD_gauche):
                    event = 'ExSkipping '
                else:
                    event = 'Backsplice '

                Val_2 = Matrice2[Tested_Junc][Case]
                q25_2 = Matrice2[Tested_Junc].quantile([0.25]).values[0]
                q50_2 = Matrice2[Tested_Junc].quantile([0.50]).values[0]
                q75_2 = Matrice2[Tested_Junc].quantile([0.75]).values[0]

                Reads_2 = Matrice2['Tot_Reads'][Case]

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

                        print(Case[0:4])
                        print(Line)
                        if nbAno == False :
                            AnoSpliceList.append(Case[0:4])
                            nbAno = True
                        AnoSpliceList.append(Line)

            l += 1

        if nbAno == True:
            AnoSpliceList.append('')

    os.chdir(rep)
    AnoSpliceList = pd.DataFrame(AnoSpliceList)
    AnoSpliceList.to_csv('Ano_' + i + '.csv', sep = ';', header = False, index = False)
