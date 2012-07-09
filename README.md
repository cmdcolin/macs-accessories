macs-accessories
================

The package Model-based analysis for ChIP-seq (MACS) is very useful software for peak finding on ChIP-seq data. The outputted files include peaks in BED and XLS formats and the 'tag pileup' in wiggle file format. We use macs-accessories for additional data analysis from the files including visualizations and calculations of normalized difference (NormDiff) scores from Zheng et al. (2010)




We created a 'wiggle class', an S3 R object, for loading MACS output files and doing additional data analysis. The class automatically loads peak files (bed,xls,csv) and and wiggle files (macs output format)



```r
wig1 = WiggleClass("S96")
wig2 = WiggleClass("HS959")
wig1$loadWiggles()
wig2$loadWiggles()
```





The normalized difference score gives us on average the expected value of the ChIP-seq subtracted from the input data using a simple random model. For $A$, $B$ representing chip-seq and input control data respectively

$A\sim Poisson(f+g)$
$B\sim Poisson(cg)$

Then the NormDiff score $Z$ is defined as

$$Z(x_i)=\frac{A(x_i)-B(x_i)/c}{\hat\sigma}$$

We use the data to estimate scaling factor $c$ and variance $\hat\sigma$

We can look at the average normalized difference scores of the peaks, and see how this compares with the same location in other experiments, simply by comparing wig1 with wig2 Z scores for the peak.




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
d <- density(datasort, adjust = 1.4)  # returns the density data
plot(d, main = "Kernel density of NormDiff scores")  # plots the results
polygon(d, col = "#BB2222CC", border = "#222244")
```

![plot of chunk zall](http://i.imgur.com/l69QG.png) 

```r
clone = datasort
qqnorm(clone)
qqline(clone, col = 2)
```

![plot of chunk zall](http://i.imgur.com/gnTF7.png) 


By taking the pvalues of all Zscores we can see which ones are significant according to a p-value $P(\bar x \leq X)$


```r
plotZscoreColor("HS959 pvalue for S96 peaks", wig1, wig2, wz1, wz2)
```

![plot of chunk color](http://i.imgur.com/UpUZ9.png) 


We want to use hypothesis testing to observe NormDiff scores that are highly different from the background. We used a likelihood ratio test sgnificance level of 1%. When comparing ChIP-seq experiments from different strains of yeast.




```r
ret = plotZscoreCutoff("S96 peaks vs HS959 synteny (Hypothesis testing)", 
    wig1, wig2, wz1, wz2, 0.05)
```

![plot of chunk cutoff](http://i.imgur.com/q8tFB.png) 

