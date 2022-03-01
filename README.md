# RT-MLPA scripts

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
