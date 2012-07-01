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





The normalized difference score gives us on average the expected value of the ChIP-seq subtracted from the input data using a simple random model. For $A$, $B$ representing chip-seq and input control data respectively

$A\sim Poisson(f+g)$
$B\sim Poisson(cg)$

Then the NormDiff score is defined as

$$Z=(A-B/c)/\hat\sigma$$

We use the data to estimate scaling factor $c$ and variance $\hat\sigma$ and then we can look at the distribution of average normalized difference scores for all peaks, and see how they correspond to syntenic regions in other datasets




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

![plot of chunk d2](figure/d21.png) 

```r
r2 = plotMaxAvgZscore("Max Avg HS959 peak NormDiff score vs S96 synteny w=100", 
    wig2, wig1, wz4, wz3, "#bb0000", "#009900")
```

![plot of chunk d2](figure/d22.png) 


If we sort the data we can get an idea of the connections between the cutoff applied to both datasets



```r
plotSortedMaxAvgZscore("Sorted HS959 Max Avg Normdiff connected with S96 peak regions", 
    wig1, wig2, r1, "#00119919", "#bb000019", "#001199", "#aa0000")
```

![plot of chunk sorted](http://i.imgur.com/ctuBS.png) 

```r
plotSortedMaxAvgZscore("Sorted S96 Max Avg Normdiff connected with HS959 peak regions", 
    wig2, wig1, r2, "#00990019", "#bb000019", "#009900", "#aa0000")
```

![plot of chunk sorted](http://i.imgur.com/ttfhR.png) 


We can use a conditional probability based on a correlation between the datasets. This is defined as, for normalized difference scores for peaks

$$ p(x|y)={{P(\bar x \leq X) corr(x,y)}\over{P(\bar y \leq Y)}} $$


If we look at colored plots of the probability we get from these plots we see something like this



```r

plotMaxAvgZscoreColor("S96 peaks vs HS959 synteny", wig1, wig2, wz1, 
    wz2)
```

![plot of chunk rainbow](http://i.imgur.com/tnBP6.png) 

```r
plotMaxAvgZscoreColor("HS959 vs S96 peaks colored by probability", 
    wig2, wig1, wz4, wz3)
```

![plot of chunk rainbow](http://i.imgur.com/GR7SR.png) 

