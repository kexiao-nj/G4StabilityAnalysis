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
RScript analysis/fig2-bcd-distr-acrs-chrmstate.R
``` 
OR you can run the R scripts line by line.

### 3.2 Causal Bayesian network Construction 1
For Stability, eG4Sig, phyloP score, ChromState, ATACSig and #TF
#### 3.2.1 Construnction of causal Bayesian networks with different sampling strategies and sizes
* Try different combinations of sampling strategies and sizes, this will take a long time (more than 24 hours) for the construction of these networks.
* **OR you can skip this step and use the pre-generated networks in the `robusttest` folder**:
```shell
RScript analysis/fig3-0-causal-test.R
```
#### 3.2.2 Robustness Measurement
* Plot coverage against accuracy:
```shell
RScript analysis/fig3-a-causal-robust-measure.R
```
#### 3.2.3 Drawing the most robust network
```shell
RScript analysis/fig3-b-causal-network-drawing.R
```


### 3.3 Causal Bayesian network Construction 2
For Stability, phyloP score, ChromState, ATACSig and Occupancy of specific TFs
#### 3.3.1 Construnction of causal Bayesian networks for different TFs
* The sampling strategy and size for the most robust network were employed, based on the test of the last step. This may take a long time for the construction of these networks.
* **OR you can skip this step and use the pre-generated networks in the `tfrobustnet` folder**:
```shell
RScript analysis/fig4-a-causal-network-constuct-tf-specific.R
```
#### 3.3.2 Drawing the networks
```shell
RScript analysis/fig4-a-causal-network-drawing-tf-specific.R
```

### 3.4 Gene Function Analysis
Biological process enrichment, with gene sets grouped by the stability of upstream G4s and cell lines:
```shell
RScript analysis/fig5-a-go-enrichment-analysis.R
```