# G4StabilityAnalysis
## 1. Description
The workflow is for the analysis of G4 stability and other related signals.

## 2. Suggested procedure for data preparation
You can prepare the data by yourself, or use the prepared datasets and skip to the causality inference step. The prepared datasets for K562, HepG2 and HEK293T cell lines are in the folder named prepared_datasets.

## 3. Analysis Procedure
### 3.0 Preparation for the output folder
```shell
if [ ! -d output-fig ]; then
    mkdir output-fig
fi
``` 
### 3.1 Multi-omics correlation analysis
For Spearman correlation calculation:
```shell
RScript analysis/fig2-a-corr.R
``` 
For the distribution of Stability, eG4Sig, ATACSig, phyloP score, and TF number across the Chromatin states:
```shell
RScript analysis/fig2-a-distr-acrs-chrmstate.R
``` 
OR you can run the R scripts line by line.

### 3.2
Will be uploaded soon.

