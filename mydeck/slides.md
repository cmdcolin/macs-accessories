---
title       : Comparison and analysis of ChIP-seq experiments
subtitle    : 
author      : Colin Diesh
job         : 
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : [mathjax]            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
---

## Outline
# Overview of the presentation

- Background
  - What is ChIP-seq?
  - Overview of problem
- Comparing ChIP-seq data
  - Review of existing methods
  - My Approach
    - Finding differential binding sites
- Results
  - Evaluate method
  - Compare with other methods
- Conclusion


--- 

## Background

# What is ChIP-seq?


- ChIP-seq = High-throughput sequencing + chromatin immunoprecipitation
- Produces millions of short DNA reads that are enriched for DNA-protein interactions


![chipseq-overview](assets/img/image1.png)
Huber et al, 2006

---
## Background

# What does ChIP-seq look like?

- It's noisey!
- Contains peaks that represent binding sites
- Noise depends on the sequencing depth, relative enrichment

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 

--- 



## Background

# Analyzing the data--finding peaks in the data

- Use knowledge about data to build peak model
- Peak model uses distribution of paired-end reads to estimate the size of the DNA fragments
  - Sequencer reads are ~20-50bp
  - Fragment size is between ~100-300bp


![chipseq-fragments](assets/img/image2.png)

Wilbanks et al, 2010


---


## Overview 
# ChIP-seq applications

- ChIP-seq can be used to learn about
  - transcription factor binding sites
  - histone modifications
  - regulatory motifs
- Important to compare datasets from experiments
- Find effects of binding sites regulation on gene expression

    
![regulatory-network](assets/img/image4.gif)
From http://genome.tbdb.org/annotation/genome/tbdb/RegulatoryNetwork.html

---

## Overview
# Comparing ChIP-seq data


- Binding sites can change between experiments, between species, etc
- Comparing related genomes can provide some insight

![odom-compare](assets/img/image8.png) 


Odom, Dowell et al. 2007

---


## Goals

- Using replicates
  - High throughput experiments generate lots of data, but replication is needed
  - Many protocols pool replicate experiments together, problems?
- Directly compare data from experiments
  - Can compare overlap of peaks, but this is indirect

--- 

## Problem- How to compare different experiments?

- MACS (Model based analysis of ChIP-seq)
  - Peak finder that compares ChIP-seq to control and uses a Poisson distribution
  - Can compare different conditions by using one as a "false control"
  - $\lambda=$ expected number of tags in scanning window
  
![poisson-distribution](assets/img/image5.png)

---

## Problem- How to compare different experiments?

- MACS
- DIME (Differential identification using mixture ensemble)
  - Directly compare two experiments as log ratios between 
  - Use mixture model to find differential 'component' of distribution

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

Different views

---
## Problem- How to compare different experiments?

- MACS
- DIME 
- NormDiff
  - Analyze many experiments at once
  - Find variable binding sites
  
  

```
##     Chr Pos   S96   S96   S96 HS959  HS959 HS959
## 1 chr01  51 11.37 30.31 22.53 25.76 110.63 77.57
## 2 chr01  61 11.37 27.28 22.53 25.76 110.63 80.13
## 3 chr01 121 11.37 24.25 21.94 20.46  56.07 73.31
## 4 chr01 131 10.23 24.25 21.94 18.19  54.56 72.46
## 5 chr01 141 10.23 22.73 20.16 18.94  50.01 65.64
## 6 chr01 151 10.23 18.19 21.94 18.94  46.98 61.38
```


Table of read scores from multiple ChIP-seq experiments

---

## Problem- How to compare different experiments?

- MACS
- DIME
- NormDiff
- Our method
  - Use t-tests to find differential binding sites in ChIP-seq data


---

## Pros and Cons

- Benefits of t-test approach
  - Includes replicates
  - Evaluates significance and size of differences

- Problems
  - Requires normalization of values
  - There is a large "Background" distribution in chip-seq that is troublesome for t-tests

![t-test](assests/img/image21.png)

Equation for Standard t-test with "pooled variance"

---

---
## Solution to normalization--Scaling values

We observe a different read depths in our experiments, so without any scaling they appear to have different distributions

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


---
## Solution to normalization--Scaling values

We perform scaling to match the median across all experiments


![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

- Benefit
  - less aggressive than quantile normalization
- Consideration
  - other algorithms for normalizing based on background distribution are being developed, for example, NCIS


```
---


## Solution to t-test

- Use a "moderated t-test" to account for background distribution
- Adds a small amount of bias estimated from the mean standard deviation across entire dataset


---

## Assumption of normality

- Tested the distributions using QQ-plots
- Log ratio of S96 vs HS959 replicates vs Standard normal

![qqplot1](assets/img/testing/image15.png)


---


## Assumption of normality

- Tested the distributions using QQ-plots
- Log ratio of two S96 replicates vs Standard normal
![qqplot2](assets/img/testing/image16.png)

---


## Performing the t-test

- The moderated t-test estimates based on log-ratios according to a linear model
  - $E[y_g]=X\alpha_g$
- Then the 

![volcano](assets/img/testing/image20.png)



---

## Comparing our results to other approaches

- We perform merging for our approach from sites that are close (<20bp away)
- Look at the rate of the merging vs the p-value thresholds

![bedtools1](assets/img/testing/image28.png)

---


## Comparing our results to other approaches

- We perform merging for our approach from sites that are close (<20bp away)
- Look at the rate of the merging vs the p-value thresholds

![bedtools2](assets/img/testing/image29.png)


---


## Results from merging

- Outperforms merging compared with results from DIME
- Has results comparable to MACS which has its own merging algorithm


![limma-macs](assets/img/testing/image39.png)

---


## Compare results with gene expression data

- We can combine ChIP-seq with gene expression data to find the influence of binding sites
- Shows the binding sites influence nearby genes

![gene-expression](assets/img/testing/image42.png)

- Light color indicates correlation with random gene
- Dark color indicates correlation with nearest gene to binding site

--- 

## Future Goals

- Evaluate different merging algorithms
  - TileMap type analysis from CisGenome (Ji et al, 2008) 
  - Use moving average merging
- Do a closer look at NormDiff
- Consider modelling background distribution differently
  - omitting background distributions (Dudoit et al, 2002)
  - modelling local variations like MACS (Zhang et al. 2008)
- Look at other types of Normalizations
  - voom from Smyth et al, used in limma
  - different types of LOESS, used in DIME
  - NCIS scaling methods, modifications to MACS


--- 

## Conclusion

- Thank you for listening!
- Questions?


