macs-accessories
================

The package Model-based analysis for ChIP-seq (MACS) is very useful software for peak finding on ChIP-seq data. The outputted files include peaks in BED and XLS formats and the 'tag pileup' in wiggle file format. We use macs-accessories for additional data analysis from the files including visualizations and calculations of normalized difference (NormDiff) scores from Zheng et al. (2010)




Typical setup preamble creates a 'wiggle class', an S3 R object and loadswiggle files based on the name parameter from MACS



```r
wig1 = WiggleClass("S96")
wig2 = WiggleClass("HS959")
wig1$loadWiggles()
wig2$loadWiggles()
####
```





The normalized difference score gives us on average the expected value of the ChIP-seq subtracted from the input data using a simple random model

$A~Poisson(f+g)$
$B~Poisson(cg)$

Then the NormDiff score is defined as

$$Z=(A-B/c)/\hat\sigma$$

We use functions to estimate scaling and variance and we can look at the distribution of average normalized difference scores for all peaks, and see how they correspond to syntenic regions in other datasets




```r
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()
wz1 = wig1$Z(wig1$peaks)
wz2 = wig2$Z(wig1$peaks)
wz4 = wig2$Z(wig2$peaks)
wz3 = wig1$Z(wig2$peaks)
r1 = plotMaxAvgZscore("Max Avg  S96 peak NormDiff score vs HS959 synteny w=100", 
    wig1, wig2, wz1, wz2, "#bb0000", "#001199")
```

```
## Error: object 'xl' not found
```

```r
r2 = plotMaxAvgZscore("Max Avg HS959 peak NormDiff score vs S96 synteny w=100", 
    wig2, wig1, wz4, wz3, "#bb0000", "#009900")
```

```
## Error: object 'xl' not found
```




If we sort the data we can get an idea of the connections between the cutoff applied to both datasets



```r
plotSortedMaxAvgZscore("Sorted HS959 Max Avg Normdiff connected with S96 peak regions", 
    wig1, wig2, r1, "#00119919", "#bb000019", "#001199", "#aa0000")
```

```
## Error: object 'r1' not found
```

```r
plotSortedMaxAvgZscore("Sorted S96 Max Avg Normdiff connected with HS959 peak regions", 
    wig2, wig1, r2, "#00990019", "#bb000019", "#009900", "#aa0000")
```

```
## Error: object 'r2' not found
```




We can use a conditional probability based on a correlation between the datasets. This is defined as, for normalized difference scores for peaks

$$ p(x|y)={{normalCDF(\bar x) corr(x,y)}\over{normalCDF(\bar y)}} $$
