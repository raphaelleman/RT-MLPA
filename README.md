# RT-MLPA scripts

## Install

Donwload scrpits

```bash
git clone https://github.com/raphaelleman/RT-MLPA.git
cd ./RT-MLPA_analysis
```

## Python packages

- pandas
- matplotlib
- sklearn
- joblib

## 14_ReadFastQV10_CoCar.py

Extract sequence from Fastq FastQ Filtres

```bash
python 14_ReadFastQV10_CoCar.py -R mySequence -O /path/to/output -I /path/to/fastq
```

## 15_Count_CoCar.py

Count fusion from sequence, need the list of probes

```bash
python 15_Count_CoCar.py -R myCount -O /path/to/output -I /path/to/mySequence -P /path/to/Probes
```

## 16_Figures_CoCar.py

Create matrices and Figures from count file, need design of probes files

```bash
python 16_Figures_CoCar.py -I /path/to/myCount -P /path/to/Probes
```

## 17_Recherche_Anos_CoCar.py

Get final analysis of fusion count

```bash
python 17_Recherche_Anos_CoCar.py -R myAnalysis -O /path/to/output -I /path/to/myCount -P /path/to/Probes
```

## Statistics analyses

Under development part to compute p-value of abnormal expression

### Input

The scripts expect to have *.count file in the following format.
```bash
Fusion;Count_Full;Count_UMI;Adj_UMI
BRCA1E22G-BRCA1E24D;13;5;5
BRCA1E12G-BRCA1E3D;1;1;1
```

Only the Adj_UMI count will be used in following analyses.

### Two step analysis

If you have control samples you can perfome analysis in two step.
First, it is the model trainning step from control samples.

```bash
Rscript ./trainFusionDistribution.r -I /path/to/myControls/ -O /path/to/myModel.RData
```

Second, the trained model is used to analyze career samples.

```bash
Rscript ./appTrainedModel.r -I /path/to/mySamples/ -m path/to/myModel.RData -O /path/to/myOutput.txt
```

### One step analysis

If you have not enougth control samples, used the cross-validation in one step analysis

```bash
Rscript ./CrossValAnalyzer.r -I /path/to/mySamples/ -O /path/to/myOutput.txt
```

#### Option for One step analysis

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
