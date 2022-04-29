# RT-MLPA pipeline

## Install:

Donwload scripts:

```bash
git clone https://github.com/raphaelleman/RT-MLPA.git
cd ./RT-MLPA
```

## Python packages:
Works with Python version 3.8.10 and use these following packages:
- pandas
- matplotlib
- sklearn
- joblib
- tqdm
- contextlib

## **Pipeline:**
The pipeline will launch every following steps explained after (ReadFastQ, Count_CoCar, Figures_CoCar, Recherche_Anos_Cocar and statistics test w/ CrossValAnalyser) from fastq file. 

Fastq file were generated using bcl2fastq tool with the following command:
``` 
singularity exec /mnt/remote_references/containers/bcl2fastq/bcl2fastq_2020-06-12_07-58-43.simg bcl2fastq --runfolder-dir /mnt/tampon/NextSeqOutput/220420_NS500765_0834_AHLK7GBGXL/ --output-dir /mnt/recherche/RNA-seq/RT-MLPA_analysis/Run_Avril/FASTQ/ --sample-sheet /mnt/recherche/RNA-seq/RT-MLPA_analysis/Run_Avril/SampleSheet_run_avril.csv --no-lane-splitting 
```

Exemple of use:
```
$ python3 RT-MLPA_pipeline.py \
    
    -I RT-MLPA_analysis/Run_Avril/FASTQ/ \
    -O RT-MLPA_analysis/Run_Avril/ \
    -P RT-MLPA_scripts/Probes_Design/ \
    -G RT-MLPA_scripts/list_genes.txt \
    -R model \
    -T 40 \
    -it 20
```

### Options:

``` -I ``` : /path/to/fastq files (expect one fastq (one file) per sample ie *--no-lane-splitting option in bcl2fastq*)

``` -O ``` : /path/to/output directory

``` -P ``` : /path/to/probes design

``` -G ``` : /path/to/gene file list in txt format

``` -R ``` : /path/to/runmode (default: align,count,figure,stat,model)

Possible modes:
- **align**: run ReadFastQV10_CoCar.py and put results in /outdir/mySequence
- **count**: run Count_CoCar.py and put results in /outdir/myCount
- **figure**: run Figures_CoCar.py and put results in /outdir/myCount/*sample*/*gene*/Figures
- **stat**: run Recherche_Anos_CoCar.py and put results in /outdir/myAnalysis
- **model**: run CrossValAnalyzer.r and put results in /outdir/myStatAnalysis_nbiter

``` -T ``` : number of threads to use for paralleling (default maximum cpu threads)

``` -it ``` : /number of iteration to perform CrossValAnalyser.r (default: 10)



## Subscripts used
### **ReadFastQV10_CoCar.py**

Extract sequence from Fastq file

```bash
python ReadFastQV10_CoCar.py -R mySequence -O /path/to/output -I /path/to/fastq
```

### **Count_CoCar.py**

Count fusion from sequence, need the list of probes

```bash
python 15_Count_CoCar.py -R myCount -O /path/to/output -I /path/to/mySequence -P /path/to/Probes -G /path/to/gene_list
```

### **Figures_CoCar.py**

Create matrices and Figures from count file, need design of probes files

```bash
python 16_Figures_CoCar.py -I /path/to/myCount -P /path/to/Probes -G /path/to/gene_list
```

### **Recherche_Anos_CoCar.py**

Get final analysis of fusion count

```bash
python 17_Recherche_Anos_CoCar.py -R myAnalysis -O /path/to/output -I /path/to/myCount -P /path/to/Probes -G /path/to/gene_list
```

### **Statistics analyses**

Under development part to compute p-value of abnormal expression

#### Input

The scripts expect to have *.count file in the following format.
```bash
Fusion;Count_Full;Count_UMI;Adj_UMI
BRCA1E22G-BRCA1E24D;13;5;5
BRCA1E12G-BRCA1E3D;1;1;1
```

Only the Adj_UMI count will be used in following analyses.

#### Two step analysis

If you have control samples you can perfom analysis in two step.
First, it is the model trainning step from control samples.

```bash
Rscript ./trainFusionDistribution.r -I /path/to/myControls/ -O /path/to/myModel.RData
```

Second, the trained model is used to analyze career samples.

```bash
Rscript ./appTrainedModel.r -I /path/to/mySamples/ -m path/to/myModel.RData -O /path/to/myOutput.txt
```

#### One step analysis

If you have not enougth control samples, used the cross-validation in one step analysis

```bash
Rscript ./CrossValAnalyzer.r -I /path/to/mySamples/ -O /path/to/myOutput.txt
```

##### Option for One step analysis

**-s, --nSamp** Integer
* Number of ramdom samples used to perform cross validation [Default=10]

**-i, --nIter** Integer
* Number of iteration perfom during cross validation [Default=10]

**-p, --plot** /path/to/plot.pdf
* Print plot of count distribution

## Authors

* Philippe Ruminy
* Raphael Leman - [raphaelleman](https://github.com/raphaelleman/ "tittle")
    * You can contact me at: r.leman@baclesse.unicancer.fr or raphael.leman@orange.fr

## License

This project is licensed under the Apache License - see the [LICENSE](https://github.com/raphaelleman/RT-MLPA/blob/master/LICENSE "tittle") file for details
