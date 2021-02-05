# RT-MLPA analysis

---

Statistical analysis of fusion count from RT-MLPA sequencing

**Table of contents**

* [Install](#1)
* [Input](#2)
* [Two step analysis](#3)
* [One step analysis](#4)
    * [Option for One step analysis](#5)
* [Authors](#6)
* [License](#7)

---

## Install<a id="1"></a>

Donwload scrpits

```bash
git clone
cd ./RT-MLPA_analysis
```

## Input<a id="2"></a>

The scripts expect to have *.count file in the following format.
```bash
Fusion;Count_Full;Count_UMI;Adj_UMI
BRCA1E22G-BRCA1E24D;13;5;5
BRCA1E12G-BRCA1E3D;1;1;1
```

Only the Adj_UMI count will be used in following analyses.

## Two step analysis<a id="3"></a>

If you have control samples you can perfome analysis in two step.
First, it is the model trainning step from control samples.

```bash
Rscript ./trainFusionDistribution.r -I /path/to/myControls/ -O /path/to/myModel.RData
```

Second, the trained model is used to analyze career samples.

```bash
Rscript ./appTrainedModel.r -I /path/to/mySamples/ -m path/to/myModel.RData -O /path/to/myOutput.txt
```

## One step analysis<a id="4"></a>

If you have not enougth control samples, used the cross-validation in one step analysis

```bash
Rscript ./CrossValAnalyzer.r -I /path/to/mySamples/ -O /path/to/myOutput.txt
```

### Option for One step analysis <a id="5"></a>

**-s, --nSamp** Integer
* Number of ramdom samples used to perform cross validation [Default=10]

**-i, --nIter** Integer
* Number of iteration perfom during cross validation [Default=10]

**-p, --plot** /path/to/plot.pdf
* Print plot of count distribution

## Authors <a id="6"></a>

* Raphael Leman - [raphaelleman](https://github.com/raphaelleman/ "tittle")
    * You can contact me at: r.leman@baclesse.unicancer.fr or raphael.leman@orange.fr

## License <a id="7"></a>

This project is licensed under the MIT License - see the [LICENSE](https://github.com/raphaelleman/SpliceLauncher/blob/master/LICENSE "tittle") file for details
