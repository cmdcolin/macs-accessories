hypothesis-testing
========================================================

Here we are working on hypothesis testing for Normalized Difference scores
In order to visualize deviations from Normal, we can use Q-Q plots and kernel density estimators






```r
wig1 = WiggleClass("S96")
wig2 = WiggleClass("HS959")
wig1$loadWiggles()
wig2$loadWiggles()
###
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()


############ Z score peaks
wz1 = wig1$Z(wig1$peaks)
wz2 = wig2$Z(wig1$peaks)
wz4 = wig2$Z(wig2$peaks)
wz3 = wig1$Z(wig2$peaks)

######### Z score all
wza1 = wig1$Zall()
wza2 = wig2$Zall()
maxw1 <- wig1$getMaxAvgZscore(wz1)
maxw2 <- wig2$getMaxAvgZscore(wz2)
wzamax1 = wig1$getMaxAvgZscoreAll(wza1)
wzamax2 = wig2$getMaxAvgZscoreAll(wza2)
```






```r

## Setup data
datasort = as.numeric(wza1[[3]][, 4])


par(mfrow = c(2, 1))



## Q-Q plot peak scores
select = (wig1$peaks[, 1] == "chr01.fsa")
datasel1 = unlist(wz1[select])
datasel2 = unlist(wz2[select])
clone = datasort
qqnorm(clone)
points(qqnorm(datasel1, plot = FALSE), col = "blue")
points(qqnorm(datasel2, plot = FALSE), col = "red")
qqline(rnorm(5000), col = 2)
legend("topleft", legend = c("Whole genome", "S96 peaks", "HS959 peaks"), 
    fill = c("black", "blue", "red"))





## Kernel density plot
d <- density(datasort, adjust = 1.4)  # returns the density data
plot(d, main = "Kernel density of NormDiff scores")  # plots the results
polygon(d, col = "#2222BBCC", border = "#222244")
d1 <- density(datasel1, adjust = 1.4)  # returns the density data
d2 <- density(datasel2, adjust = 1.4)  # returns the density data
lines(d1)
lines(d2)
polygon(d1, col = "#BB222288", border = "#222244")
polygon(d2, col = "#22BB2288", border = "#222244")
legend("topright", legend = c("Whole genome", "S96 peaks", "HS959 peaks"), 
    fill = c("blue", "green", "red"))
```

![plot of chunk test](figure/test1.png) 

```r



par(mfrow = c(2, 1))

## Q-Q plot MEAN peak scores
dataselmean = wzamax1[, 4]
dataselmean1 = maxw1
dataselmean2 = maxw2
clone = dataselmean
qqnorm(clone)
points(qqnorm(dataselmean1, plot = FALSE), col = "blue")
points(qqnorm(dataselmean2, plot = FALSE), col = "red")
qqline(clone, col = 2)
legend("topleft", legend = c("Whole genome", "S96 peaks", "HS959 peaks"), 
    fill = c("black", "blue", "red"))





## Kernel density plot MEANS
d <- density(dataselmean, adjust = 1.4)  # returns the density data
```

```
## Error: 'x' contains missing values
```

```r
plot(d, main = "Kernel density of MEAN NormDiff scores")  # plots the results
polygon(d, col = "#2222BBCC", border = "#222244")
d1 <- density(dataselmean1, adjust = 1.4)  # returns the density data
d2 <- density(dataselmean2, adjust = 1.4)  # returns the density data
lines(d1)
lines(d2)
polygon(d1, col = "#BB222288", border = "#222244")
polygon(d2, col = "#22BB2288", border = "#222244")
legend("topright", legend = c("Whole genome", "S96 peaks", "HS959 peaks"), 
    fill = c("blue", "green", "red"))
```

![plot of chunk test](figure/test2.png) 

```r


```



       
