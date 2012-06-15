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


plot(r1[['2shared']],dnorm(r1[['2shared']]),pch='.',col='blue',xlim=c(-10,10),ylim=c(0,0.4))
points(r1[['2unique']],dnorm(r1[['2unique']]),pch='.',col='red')

plot(r1[['2shared']],dexp(r1[['2shared']]),pch='.',col='blue',xlim=c(-10,10),ylim=c(0,0.4))
points(r1[['2unique']],dexp(r1[['2unique']]),pch='.',col='red')

plot(r1[['2shared']],dt(r1[['2shared']],1),pch=1,col='blue',xlim=c(-10,10),ylim=c(0,0.4))
points(r1[['2unique']],dt(r1[['2unique']],1),pch='.',col='red')
plot(r1[['2shared']],pt(r1[['2shared']],10),pch='.',col='blue')
points(r1[['2unique']],pt(r1[['2unique']],10),pch='.',col='red')
       
       
       

r2s=c(r1[['2shared']],r1[['2unique']])
corr1=match(r1[['2shared']],r2s)
corr2=match(r1[['2unique']],r2s)
r2s=(r2s-mean(r2s))/sqrt(var(r2s))
r2snorm=dnorm(r2s,mean(r2s),sqrt(var(r2s)))
r2s1norm=dnorm(r1[['2shared']],mean(r2s),sqrt(var(r2s)))
r2s2norm=dnorm(r1[['2unique']],mean(r2s),sqrt(var(r2s)))
plot(r2s,r2snorm,pch='.',col='blue')
points(r2s[corr1],r2snorm[corr1],pch=20,col='red')
points(r2s[corr2],r2snorm[corr2],pch=20,col='green')


library(xtable)
debug=TRUE
env=loadMacsEnv('S96rep1','S96rep2')
attach(env)
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()
xtable(cbind(c(wig1$scaling,wig1$variance),c(wig2$scaling,wig2$variance))) 
wig1$peaks
wig1$shared
wig1$unique
