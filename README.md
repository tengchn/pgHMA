# pgHMA: Application of the Heteroduplex Mobility Assay Analysis in Phylogenetics and Population Genetics.

This pipeline is used to perform the phylogenetics and population genetic analyses based on heteroduplex mobility data, with the comparison of related analyses from DNA sequences.

**GetPeaks**: automatic calculate the heteroduplex mobility distance from the results of LabChip GX software. *Contributed by Thomas Wong.*

**pgHMAtools**: R package including some useful tools for pgHMA analyses.

**HMA_pairwise_analysis.R**: R script for phylogenetic analysis with bootstrapping and skyline plots to estimate the demographic history of a population using heteroduplex mobility distances, as well as comparing these results with the same analyses based on nucleotide sequences.

**HMA_mismatch_optimization.R**: R script for optimizing frequency distributions of heteroduplex mobility distances to minimize the differences of mismatch distributions obtained using nucleotide sequences, and plot related figures for each mixture with before and after optimization.

## Requirements
- g++
- RStudio

## Install
```
git clone https://github.com/tengchn/pgHMA.git
cd GetPeaks
#To compile the program of getPeaks
g++ -o getPeaks getPeaks.cpp
```
```
#To install the R package of pgHMAtools in R
devtools::install_github("tengchn/pgHMA/pgHMAtools")
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
**HMA_pairwise_analysis.R and HMA_mismatch_optimization.R**
```
Usage: run the scripts for all analyses in RStudio.
```
