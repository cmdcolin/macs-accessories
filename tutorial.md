macs-accessories
================

The package Model-based analysis for ChIP-seq (MACS) is very useful software for peak finding on ChIP-seq data. The outputted files include peaks in BED and XLS formats and the 'tag pileup' in wiggle file format. We use macs-accessories for additional data analysis from the files including visualizations and calculations of normalized difference (NormDiff) scores from Zheng et al. (2010)




Typical setup preamble creates a 'wiggle class', an S3 R object and loadswiggle files based on the name parameter from MACS



```r
wig1 = WiggleClass("S96")
wig2 = WiggleClass("HS959")
wig1$loadWiggles(globalenv())
wig2$loadWiggles(globalenv())
####
```





The normalized difference score gives us on average the expected value of the ChIP-seq subtracted from the input data using a simple random model

$A~Poisson(f+g)$
$B~Poisson(cg)$

Then the NormDiff score is defined as

$Z=(A-B/c)/\hat\sigma$

We can look at the distribution of average normalized difference scores on the plots




```r
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()
wz1 = wig1$Z(wig1$peaks)
wz2 = wig2$Z(wig1$peaks)
wz4 = wig2$Z(wig2$peaks)
wz3 = wig1$Z(wig2$peaks)
r1 = plotMaxAvgZscore("Max Avg S96 peak NormDiff score vs HS959 synteny w=100", 
    wig1, wig2, wz1, wz2, "lightblue", "orange")
```

![plot of chunk d2](http://i.imgur.com/bUNTU.png) 

```r
r2 = plotMaxAvgZscore("Max Avg HS959 peak NormDiff score vs S96 synteny w=100", 
    wig2, wig1, wz4, wz3, "pink", "orange")
```

![plot of chunk d2](http://i.imgur.com/jZZaK.png) 


If we sort the data we can get an idea of the connections between the cutoff applied to both datasets



```r
plotSortedMaxAvgZscore("Sorted HS959 Max Avg Normdiff connected with S96 peak regions", 
    wig1, wig2, r1, "#00119919", "#bb000019", "#001199", "#aa0000")
```

![plot of chunk sorted](http://i.imgur.com/1GDeM.png) 

```r
plotSortedMaxAvgZscore("Sorted S96 Max Avg Normdiff connected with HS959 peak regions", 
    wig2, wig1, r2, "#00990019", "#bb000019", "#009900", "#aa0000")
```

![plot of chunk sorted](http://i.imgur.com/sP2pW.png) 


