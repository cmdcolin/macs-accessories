source('readplot.R')
debug=FALSE
wigenv=new.env()
wig1=WiggleClass('S96')
wig2=WiggleClass('HS959')
wig1$loadWiggles(wigenv) 
wig2$loadWiggles(wigenv)
###
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()


############
# Z score peaks
wz1=wig1$Z(wig1$peaks)
wz2=wig2$Z(wig1$peaks)
wz4=wig2$Z(wig2$peaks)
wz3=wig1$Z(wig2$peaks)




#########
# Z score all
wza1=wig1$Zall()
wza2=wig2$Zall()
maxw1<-wig1$getMaxAvgZscore(wz1)
maxw2<-wig2$getMaxAvgZscore(wz2)
wzamax1=wig1$getMaxAvgZscoreAll(wza1)
wzamax2=wig2$getMaxAvgZscoreAll(wza2)



#####################3
# SCRAP CODE

#wig1$loadWiggles() 
#wig2$loadWiggles()
###### densityplot
#library(lattice)
#d<-densityplot(~datasort,data=faithful,kernel="rect")
#d<densityplot(~datasort,kernel="gaussian")
#plot(d, main='Kernel density of NormDiff scores') # plots the results

##############
# Plot Overlaps
#bedselect=wig2$peaks[b1,]
#selection1=wig1$Z(bedselect)
#selection2=wig2$Z(bedselect)
##reads1=wig1$getChipReads(bedselect)
#reads2=wig2$getChipReads(bedselect)


#### Draw rainbow
#plot(c(100, 250), c(300, 450), type = "n",main="r")
#i <- 4*(0:10)
## draw rectangles with bottom left (120, 300)+i  and top right (177, 380)+i
#rect(120, 300+i, 177, 380+i, col=rainbow(11, start=.0,end=.7))

#plot(p1,type='l',col="#aa0000bb",ylim=c(0,20),xlim=c(0,100000))
#lines(p2,col="#101010cc")



#plotZall2(c1,wig1,wig2)
#plotZall3(c1,wig1,wig2)

#plotMaxAvgReads('S96 avg peak reads vs HS959 syntenic', wig1, wig2, 'green', 'blue') 
#plotMaxAvgReads('HS959 avg peak reads vs S96 syntenic', wig2, wig1, 'green', 'red')
