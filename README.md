# pgHMA: Application of the Heteroduplex Mobility Assay Analysis in Phylogenetics and Population Genetics.

This pipeline is used to perform the phylogenetics and population genetic analyses based on heteroduplex mobility data, with the comparison of related analyses from DNA sequences.

**GetPeaks**: automatic calculate the heteroduplex mobility distance from the results of LabChip GX software.

**HMA_pairwise_analysis.R**:  

## Requirements
- g++
- RStudio 
- Geneious 
- BEAST 
- Tracer 
- IQ-TREE

## Install
```
git clone https://github.com/tengchn/pgHMA.git
cd GetPeaks
#To compile the program
g++ -o getPeaks getPeaks.cpp
```

## Usage
**GetPeaks**
```
[Progarm for automatic calculating the heteroduplex mobility distance from the results of LabChip GX software]

To get the peaks from the input data:
  Syntax: ./getPeaks view [csv file] (options)

  Options:
     -k [integer value of k, default: 0]
     -l [integer value of l, default: 1]
     -t [double value of t, default: 0.6]

  Usage: ./getPeaks view Data/Raw data from LabChip GXII/LabChip output for 16Kangaroos pairwise combination.csv

To identify the homo-duplex and hetero-duplex, assuming it is a mixture of TWO samples:
  Syntax: ./getPeaks compute [csv file] -m [time (in double) of lower marker] -n [time (in double) of upper marker] (options)

  Options:
     -k [integer value of k, default: 0]
     -l [integer value of l, default: 1]
     -tlist [list of t values, default: 0.3,0.6,1.0]

  Usage: ./getPeaks compute Data/Raw data from LabChip GXII/LabChip output for 16Kangaroos pairwise combination.csv -m 11.92 -n 22.65

NOTE: If some samples have errors in the output file, please visualize and check these samples with LabChip GX software as there may have erroneous peaks caused by signal noise.
```
