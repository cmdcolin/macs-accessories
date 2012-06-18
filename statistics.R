library(xtable)
debug=TRUE
env=loadMacsEnv('S96','HS959')
attach(env)
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()
xtable(cbind(c(wig1$scaling,wig1$variance),c(wig2$scaling,wig2$variance))) 



wz1=wig1$Z(wig1$peaks)
wz2=wig2$Z(wig1$peaks)
wz4=wig2$Z(wig2$peaks)
wz3=wig1$Z(wig2$peaks)

r1=plotMaxAvgZscore('Mean S96 peak Normdiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2,'yellow', 'darkviolet')
r2=plotMaxAvgZscore('Mean HS959 peak Normdiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'orange','darkviolet')
r3=plotAvgZscore('Mean S96 peak NormDiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2, 'lightblue', 'orange')
r4=plotAvgZscore('Mean HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'pink','orange')

plotSortedMaxAvgZscore('Sorted HS959 Max Avg Normdiff in S96 peak regions',wig1,wig2,r1,r2,'red','blue')
plotSortedMaxAvgZscore('Sorted S96 Normdiff in HS959 peak regions',wig2,wig1,r2,r1,'green','blue')
plotSortedMaxAvgZscore('Sorted S96 peak Avg NormDiff score vs HS959 synteny w=100', wig1, wig2, r3, r4, 'red', 'darkblue')
plotSortedMaxAvgZscore('Sorted HS959 peak Avg NormDiff score vs S96 synteny w=100', wig2, wig1, r4, r3, 'green','darkblue')



plotSortedMaxAvgZscoreX('Sorted HS959 Max Avg Normdiff in S96 peak regions',wig1,wig2,r1,'red','blue')
plotSortedMaxAvgZscoreX('Sorted S96 Normdiff in HS959 peak regions',wig2,wig1,r2,'green','blue')
plotSortedMaxAvgZscoreX('Sorted S96 peak Avg NormDiff score vs HS959 synteny w=100', wig1, wig2, r3, 'red', 'darkblue')
plotSortedMaxAvgZscoreX('Sorted HS959 peak Avg NormDiff score vs S96 synteny w=100', wig2, wig1, r4, 'green','darkblue')



library(xtable)
debug=TRUE
env=loadMacsEnv('S96','HS959')
attach(env)
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()
xtable(cbind(c(wig1$scaling,wig1$variance),c(wig2$scaling,wig2$variance))) 
wig1$peaks
wig1$shared
wig1$unique



debug=TRUE
wig1=WiggleClass('S96')
wig2=WiggleClass('HS959')
wig1$loadWiggles(globalenv()) 
wig2$loadWiggles(globalenv())
###
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()
Rprof()
wza1=wig1$Zall()
Rprof(NULL)
summaryRprof()
wz1=wig1$Z(wig1$peaks)
wz2=wig2$Z(wig1$peaks)
wz4=wig2$Z(wig2$peaks)
wz3=wig1$Z(wig2$peaks)
r1=plotMaxAvgZscore('Max Avg S96 peak NormDiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2,'lightblue', 'orange')
r2=plotMaxAvgZscore('Max Avg HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'pink','orange')
r3=plotAvgZscore('Mean S96 peak NormDiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2, 'lightblue', 'orange')
r4=plotAvgZscore('Mean HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'pink','orange')

plotSortedMaxAvgZscoreX('Sorted HS959 Avg Normdiff in S96 peak regions',wig1,wig2,r3,'red','blue')
plotSortedMaxAvgZscoreX('Sorted S96 Avg Normdiff in HS959 peak regions',wig2,wig1,r4,'green','blue')
plotSortedMaxAvgZscoreX('Sorted HS959 Max Avg Normdiff in S96 peak regions',wig1,wig2,r1,'red','blue')
plotSortedMaxAvgZscoreX('Sorted S96 Max Avg Normdiff in HS959 peak regions',wig2,wig1,r2,'green','blue') 




plotSortedMaxAvgZscore('Sorted HS959 Max Avg Normdiff connected with S96 peak regions',wig1,wig2,r1,'#00334455','#55221122')
plotSortedMaxAvgZscore('Sorted S96 Max Avg Normdiff connected with HS959 peak regions',wig2,wig1,r2,'#00661155','#55221122') 


