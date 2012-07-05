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

![plot of chunk d2](http://i.imgur.com/R39TV.png) 

```r
r2 = plotMaxAvgZscore("Max Avg HS959 peak NormDiff score vs S96 synteny w=100", 
    wig2, wig1, wz4, wz3, "#bb0000", "#009900")
```

![plot of chunk d2](http://i.imgur.com/zwUUU.png) 


We can calculate NormDiff scores for the whole genome using 


```r
wza1 = wig1$Zall()
wza2 = wig2$Zall()
```




We can observe the distribution of NormDiff scores



```r
datasort = as.numeric(wza1[[1]][, 4])
d <- density(datasort, adjust = 1.2)  # returns the density data
plot(d, main = "Kernel density of NormDiff scores")  # plots the results
polygon(d, col = "#BB2222CC", border = "#222244")
```

![plot of chunk zall](http://i.imgur.com/xJz4N.png) 

```r
clone = datasort
qqnorm(clone)
qqline(clone, col = 2)
```

![plot of chunk zall](http://i.imgur.com/KEQif.png) 


We want to use hypothesis testingto obsereve transcriptionfactor binding sites onthetail ofthe distribution. Plotting the max average normdiff scores from peak regions shows that these are considerably different from the population.



```r
wzamax1 = wig1$getMaxAvgZscoreAll(wza1)
wzamax2 = wig2$getMaxAvgZscoreAll(wza2)
datamax1 = wzamax1[, 4]
datamax2 = wzamax2[, 4]
qqnorm(datasort)
points(qqnorm(datamax1, plot = FALSE), col = "blue")
points(qqnorm(datamax2, plot = FALSE), col = "red")
qqline(rnorm(5000), col = 2)
```

![plot of chunk test](http://i.imgur.com/XnPIl.png) 

