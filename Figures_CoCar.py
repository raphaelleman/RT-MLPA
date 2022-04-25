import os
from re import S
import time
import argparse
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot
from matplotlib import patches
from matplotlib import transforms
from matplotlib.patches import Rectangle
from joblib import Parallel, delayed
pd.set_option('mode.chained_assignment', None)
import warnings
warnings.filterwarnings('ignore')

def get_arguments():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('-P', '--probes', dest='probes', type=str,
                        help='/path/to/probes ref')
    parser.add_argument('-I', '--input', dest='input', type=str,
                        help='/path/to/count directory')
    parser.add_argument('-T', '--threads', dest='threads', type=str, default=os.cpu_count(),
                        help='number of threads to use for paralleling (default maximum cpu threads)')
    parser.add_argument('-G', '--genes_file', dest='gene_file', type=str,
                        help='/path/to/genes file')
    return parser.parse_args()  


def draw_figure(sample, genes_list, counts_dir, probes_dir):
    #print(sample)
    for genes in genes_list:
        #print(genes)
        ## Splice matrix generation
        # Probes loading
        os.chdir(probes_dir)
        probe_name = genes + '_Probes.csv'
        probes = pd.read_csv(probe_name, sep=';' , encoding='latin-1')

        probes['DG'] = '_'
        probes['10bp'] = '_'
        probes['S'] = '0'
        for i in probes.index:
            probes['DG'][i] = (probes['Sonde'][i])[len(probes['Sonde'][i])-1]
            probes['10bp'][i] = (probes['Seq'][i])[0:10]
            probes['S'][i] = len(probes['Seq'][i])

        left_probes = probes[probes['DG'] == 'G']
        right_probes = probes[probes['DG'] == 'D']
        left_probes = left_probes.reset_index(drop = True)
        right_probes = right_probes.reset_index(drop = True)

        os.chdir(counts_dir + os.path.sep + sample + os.path.sep +  genes + os.path.sep + 'Counts')
        # Matrix generation
        countFile = os.listdir(os.getcwd())
        splices = pd.DataFrame(columns = right_probes['Sonde'], index = left_probes['Sonde'])
        if len(countFile) != 0:
            for file in os.listdir(os.getcwd()):
                os.chdir(counts_dir + os.path.sep + sample + os.path.sep +  genes + os.path.sep + 'Counts')
                Counts = pd.read_csv(file, sep=';' , encoding='latin-1')

                for x in right_probes['Sonde']:
                    for y in left_probes['Sonde']:
                        for i in range(0,len(Counts)):
                            if (x in Counts['Fusion'][i]) == True and (y in Counts['Fusion'][i]) == True:
                # Correction par UMI.... ou pas?
                                #splices[x][y] = Counts['Adj_UMI'][i]
                #splices = splices.fillna(0)
                                splices[x][y] = Counts['Count_Full'][i]
                splices = splices.fillna(0)
                # Save
                os.chdir(counts_dir+ os.path.sep + sample + os.path.sep +  genes + os.path.sep + 'Figures/Matrices')
                splices.to_csv(file, sep = ';', header = True)

        else:
            for x in right_probes['Sonde']:
                for y in left_probes['Sonde']:
                    splices[x][y] = 0
        ### Figures
        os.chdir(probes_dir)
        gene = pd.read_csv(genes + '_Design.csv', sep=';' , encoding='latin-1')

        # Graph _ mise echelle automatique fonction fichier design
        #print("graph")
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

        # Gene structure graph
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

        # Splices graphs
        spliceval = np.array(splices.loc[:,list(splices)[1]:])
        lower = np.min(spliceval)
        higher = np.max(spliceval) * 1.1
        if higher == 0:
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

        total_reads = sum(Counts['Count_Full'])
        if sample != '#Moy':
            total_UMI = sum(Counts['Adj_UMI'])
            if total_UMI != 0:
                AmpF = round(total_reads/total_UMI*10)/10
            else:
                AmpF = 0
            total_reads = sample + '  ' + genes + '     Nbr Reads : ' + str(total_reads) + '  Nbr Lig : ' + str(total_UMI) + '   Ampli : ' + str(AmpF)
        else:
            total_reads = sample + '  ' + genes + '     Nbr Reads : ' + str(total_reads)
        axes.text(50, 1,total_reads, fontsize = 8)

        filename =  sample + '_' + '_Fig_CFull_' + genes + '.png'
        os.chdir(counts_dir + os.path.sep + sample + os.path.sep +  genes + os.path.sep + 'Figures')
        figure.savefig(filename, bbox_inches = 'tight', dpi = 400)
        pyplot.close('all')
        figure.clear('all')



def main():
    start_time = time.time()

    #get arguments
    args = get_arguments()
    probes_dir = os.path.realpath(args.probes)
    counts_dir = os.path.realpath(args.input)
    gene_file = os.path.realpath(args.gene_file)
    threads_number = args.threads
    
    #get sample list
    samples_list = [file for file in os.listdir(counts_dir)]

    # get gene list from file
    genes_list = []
    with open(gene_file, "r") as file:
        for line in file:
            genes_list.append(line.strip())
    print(genes_list)

    #launch analysis
    Parallel(n_jobs= threads_number, verbose = 1)(delayed(draw_figure)(sample, genes_list, counts_dir, probes_dir) for sample in samples_list)

    print(f"--- Program executed in {float(time.time() - start_time) / 60} minutes ---")



if __name__ == "__main__":
    main()
