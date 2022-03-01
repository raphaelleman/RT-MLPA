import os
import statistics
import shutil
import csv
import math
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot
from matplotlib import patches
from matplotlib import transforms
from matplotlib.patches import Rectangle
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import argparse
pd.set_option('mode.chained_assignment', None)

#run = '2021_04_14_CoCar1'
#run = '2021_06_04_CoCar2'
#run = '2021_09_15_CoCar3'
my_parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

my_parser.add_argument('-P','--probes',type=str,
                       help='/path/to/probes ref')

my_parser.add_argument('-I','--input',type=str,
                       help='/path/to/count directory')

args = my_parser.parse_args()

#### Defnition Analyses
Analysis = 'BRCA'
#AnalysisSp = 'BRCA1'

#run = '2021_04_14_CoCar1'
#run = '2021_06_04_CoCar2'
#run = '2021_09_15_CoCar3'
#run = '2021_09_27_CoCar4'
#run = 'Ok_Cocar'
repProbes = os.path.realpath(args.probes)

### Listes FastQ et Echantillons
repCount = os.path.realpath(args.input)
Sample = []

for file in os.listdir(repCount):
    if file[0:7] != 'Matrice' and file[0:3] != 'Ano':
        Sample.append(file)

#CoCar_Genes = ['APC','AXIN2','BMPR1A','BUB1','EPCAM','FAN1','GALNT12','GREM1','MLH1','MSH2','MSH3','MSH6','MUTYH','NTHL1','PMS2','POLD1','POLE','PTEN','RNF43','RPS20','SMAD4','STK11','TP53']
#CoCar_Genes = ['APC','AXIN2','BMPR1A','BUB1','EPCAM','FAN1','GALNT12','MLH1','MSH2','MSH3','MSH6','MUTYH','NTHL1','PMS2','POLD1','POLE','PTEN','RNF43','RPS20','SMAD4','STK11','TP53']
#CoCar_Genes = ['APC','BMPR1A','EPCAM','MLH1','MSH2','MSH6','MUTYH','NTHL1','PMS2','POLD1','POLE','PTEN','SMAD4','STK11','TP53','CDH1','PALB2','RAD51C','BRCA1','BRCA2','RAD51D']

CoCar_Genes = ['BRCA1','BRCA2','PALB2','RAD51C','RAD51D']

##### Analyse
for sp in Sample:
    print(sp)
    for Analysis in CoCar_Genes:
        print(Analysis)
        ## Generation des Matrices de splice
        # Chargement des sondes
        os.chdir(repProbes)
        probes = Analysis + '_Probes.csv'
        sondes = pd.read_csv(probes, sep=';' , encoding='latin-1')

        sondes['DG'] = '_'
        sondes['10bp'] = '_'
        sondes['S'] = '0'
        for i in sondes.index:
            sondes['DG'][i] = (sondes['Sonde'][i])[len(sondes['Sonde'][i])-1]
            sondes['10bp'][i] = (sondes['Seq'][i])[0:10]
            sondes['S'][i] = len(sondes['Seq'][i])

        sondesG = sondes[sondes['DG'] == 'G']
        sondesD = sondes[sondes['DG'] == 'D']
        sondesG = sondesG.reset_index(drop = True)
        sondesD = sondesD.reset_index(drop = True)

        os.chdir(repCount + os.path.sep + sp + os.path.sep +  Analysis + os.path.sep + 'Counts')
        # Generation de la Matrice
        countFile = os.listdir(os.getcwd())
        splices = pd.DataFrame(columns = sondesD['Sonde'], index = sondesG['Sonde'])
        if len(countFile) != 0:
            for file in os.listdir(os.getcwd()):
                os.chdir(repCount + os.path.sep + sp + os.path.sep +  Analysis + os.path.sep + 'Counts')
                Counts = pd.read_csv(file, sep=';' , encoding='latin-1')

                for x in sondesD['Sonde']:
                    for y in sondesG['Sonde']:
                        for i in range(0,len(Counts)):
                            if (x in Counts['Fusion'][i]) == True and (y in Counts['Fusion'][i]) == True:
                # Correction par UMI.... ou pas?
                                #splices[x][y] = Counts['Adj_UMI'][i]
                #splices = splices.fillna(0)
                                splices[x][y] = Counts['Count_Full'][i]
                splices = splices.fillna(0)
                # Sauvegarde
                os.chdir(repCount+ os.path.sep + sp + os.path.sep +  Analysis + os.path.sep + 'Figures/Matrices')
                splices.to_csv(file, sep = ';', header = True)

        else:
            for x in sondesD['Sonde']:
                for y in sondesG['Sonde']:
                    splices[x][y] = 0
        ### Figures
        os.chdir(repProbes)
        gene = pd.read_csv(Analysis + '_Design.csv', sep=';' , encoding='latin-1')

        # Graphe _ mise echelle automatique fonction fichier design
        figure = pyplot.figure()
        axes = figure.add_subplot(111)

        k = (gene.loc[3]).tolist()
        k[0] = 0
        k = [int(i) for i in k]
        k = max(k)+ 10
        axes.set_xlim(0, k)

        k = (gene.loc[4]).tolist()
        k[0] = 0
        k = [int(i) for i in k]
        k = max(k)+ 50
        axes.set_ylim(0, k)
        axes.set_frame_on(False)
        axes.xaxis.set_visible(False)
        axes.yaxis.set_visible(False)

        # Graphique Structure des Genes
        for i in range(1, gene.shape[1]):
            Num = str(gene.iloc[0,i])
            Type = gene.iloc[1,i]
            Size = int(gene.iloc[2,i])
            XPos = int(gene.iloc[3,i])
            Ypos = 60

            if Type == 'I':
                axes.add_artist(matplotlib.lines.Line2D((XPos, XPos+Size), (Ypos, Ypos),linewidth = 0.5, color = 'black'))
            if Type == 'G':
                axes.add_artist(patches.Rectangle((XPos,Ypos-2), Size, 4,edgecolor = 'black', facecolor = 'black',fill = False, linestyle = 'solid',linewidth = 0.5, zorder = 1))
            if Type == 'D' and Size != 0 :
                axes.add_artist(patches.Rectangle((XPos,Ypos-2), Size, 4,edgecolor = 'black', facecolor = 'black',fill = False, linestyle = 'solid',linewidth = 0.5, zorder = 1))

            Fz = 4
            Cor = 1

            if Num != '0'and len(Num) < 3:
                axes.text(XPos+Size/2,Ypos - 5, Num, fontsize = Fz)

            if Num != '0'and len(Num) > 3:
                axes.text(XPos+Size/2,Ypos - Cor, Num, rotation = 90, fontsize = Fz)

        # Graphique Splices
        spliceval = np.array(splices.loc[:,list(splices)[1]:])
        lower = np.min(spliceval)
        higher = np.max(spliceval) * 1.1
        if higher==0:
            higher = 1
        for droite in splices.columns:
                for gauche in splices.index:
                        YSplice = splices[droite][gauche]
                        start = int(gene[gauche][3]) + int(gene[gauche][2])
                        stop = int(gene[droite][3]) + int(gene[droite][2])
                        YposStart = int(gene[gauche][4])
                        YposStop = int(gene[droite][4])
                        Nat = gene[gauche][1]
                        NatD = gene[droite][1]

                        # Usual Splices
                        if stop >= start and YposStart != 0 and Nat != 'S':
                            coul = 'green'
                            ep = 0.4
                            if Nat == 'Geno':
                                coul = 'blue'
                                ep = 0.5
                            if YSplice >= (higher * 0.2 / 100):
                                high = YposStart+2 + (YSplice / higher * 40)
                                axes.add_artist(matplotlib.lines.Line2D((start, start+(stop-start)/2), (YposStart+2, high), color = coul,linewidth = ep))
                                axes.add_artist(matplotlib.lines.Line2D((start+(stop-start)/2, stop), (high, YposStart+2), color = coul,linewidth = ep))

                        # SNPs
                        if stop >= start and YposStart != 0 and Nat == 'S':
                            coul = 'red'
                            ep = 1.3
                            if NatD == 'N':
                                coul = 'green'
                            if YSplice >= (higher * 0.5 / 100):
                                high = YposStart+2 + (YSplice / higher * 45)
                                axes.add_artist(matplotlib.lines.Line2D((stop, stop), (YposStart+2, high), color = coul,linewidth = ep))

                        # Back Splices
                        if start >= stop :
                            ep = 0.3
                            if YSplice >= (higher * 0.5 / 100):
                                high = YposStart-2 - (YSplice / higher * 40)
                                axes.add_artist(matplotlib.lines.Line2D((start, start+(stop-start)/2), (YposStart-2, high), color = 'red',linewidth = ep))
                                axes.add_artist(matplotlib.lines.Line2D((start+(stop-start)/2, stop), (high, YposStop-2), color = 'red',linewidth = ep))

        TotReads = sum(Counts['Count_Full'])
        if sp != '#Moy':
            TotUMI = sum(Counts['Adj_UMI'])
            if TotUMI != 0:
                AmpF = round(TotReads/TotUMI*10)/10
            else:
                AmpF = 0
            TotReads = sp + '  ' + Analysis + '     Nbr Reads : ' + str(TotReads) + '  Nbr Lig : ' + str(TotUMI) + '   Ampli : ' + str(AmpF)
        else:
            TotReads = sp + '  ' + Analysis + '     Nbr Reads : ' + str(TotReads)
        axes.text(50, 1,TotReads, fontsize = 8)

        nomfichier =  sp + '_' + '_Fig_CFull_' + Analysis + '.png'
        os.chdir(repCount + os.path.sep + sp + os.path.sep +  Analysis + os.path.sep + 'Figures')
        figure.savefig(nomfichier, bbox_inches = 'tight', dpi = 400)
        pyplot.close('all')
        figure.clear('all')
