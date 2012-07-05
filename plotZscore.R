

r1=plotMaxAvgZscore('Max Avg  S96 peak NormDiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2,'orange', 'blue')
r2=plotMaxAvgZscore('Max Avg HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'orange','darkviolet')



r3=plotAvgZscore('Mean S96 peak NormDiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2, 'lightblue', 'orange')
r4=plotAvgZscore('Mean HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'pink','orange')



plotSortedMaxAvgZscoreX('Sorted HS959 Avg Normdiff in S96 peak regions',wig1,wig2,r3,'red','blue')
plotSortedMaxAvgZscoreX('Sorted S96 Avg Normdiff in HS959 peak regions',wig2,wig1,r4,'green','blue')
plotSortedMaxAvgZscoreX('Sorted HS959 Max Avg Normdiff in S96 peak regions',wig1,wig2,r1,'red','blue')
plotSortedMaxAvgZscoreX('Sorted S96 Max Avg Normdiff in HS959 peak regions',wig2,wig1,r2,'green','blue') 





plotSortedMaxAvgZscore('Sorted HS959 Max Avg Normdiff connected with S96 peak regions',wig1,wig2,r1,'#00334407','#55221133') 
plotSortedMaxAvgZscore('Sorted S96 Max Avg Normdiff connected with HS959 peak regions',wig2,wig1,r2,'#00661107','#55221133')
plotSortedMaxAvgZscore('Sorted HS959 Max Avg Normdiff connected with S96 peak regions',wig1,wig2,r1,'#00334407','#55221133','#003344','#552211') 
plotSortedMaxAvgZscore('Sorted S96 Max Avg Normdiff connected with HS959 peak regions',wig2,wig1,r2,'#00661107','#55221133','#006611','#552211')





plotBarCharts(wig1,wig2,wz1,wz2)

unique=uniqueBed(wig1$peaks,wig2$peaks)
shared=intersectBed(wig1$peaks,wig2$peaks)
uniqueid=match(unique$V4,wig1$peaks$V4)
sharedid=match(shared$V4,wig1$peaks$V4)  
maxw1<-wig1$getMaxAvgZscore(wz1)
maxw2<-wig2$getMaxAvgZscore(wz2)






############ plotZ All
c1=plotZall(wzamax1,wig1,wig2)
c2=plotZall(wzamax2,wig2,wig1)
plotZall(wzamax1,wig1,wig2,c1)
plotZall(wzamax2,wig2,wig1,c2)
