################
# Main file for analyzezscore-joinmod.R 
#






## See histogram of average read depth from all experiments
dist1<-apply(wiggleTable[,c(-1,-2)],2,mean)
hist(dist1)



#resize for heatmap
resize.win(10,30)



# Conduct heatmap (omit chr and pos columns)
doheatmap(wiggleTable[,c(-1,-2)],1000)




# Background subtraction and scaling
for(i in 1:nsamples) {
  r1<-paste0('V',i+1)
  r2<-paste0('V',i+2)
  ratio=median(table[[r1]]/table[[r2]])
  if(ratio>1) {
    printf("scale down control %f\n", ratio)
    table[[r1]]<-table[[r1]]/ratio
  } else {
    printf("scale down treat %f\n", 1/ratio)
    table[[r2]]<-table[[r2]]*ratio
  }
}

tablescale<-table

for(i in 1:length(dist1)) {
  r1<-paste0('V',i)
  printf("Scale %s by mean total read depth %f\n", r1, dist1[i])
  
  tablescale[[r1]]<-table[[r1]]/dist1[i]
}






table<-wiggleTable
backSubList<-lapply(1:nsamples,function(i) {
  pos=i*2
  str1<-paste0('V',pos-1)
  str2<-paste0('V',pos)
  control<-table[[str1]]
  treat<-table[[str2]]
  
  m1=median(control/treat)
  treat-control/m1
})

# combine rows into table
backSubTable<-as.data.frame(do.call(cbind,backSubList))

#get pos and chr columns
backSubTable<-with(wiggleTable, cbind(pos,chr,backSubTable))


#Comparison plots
source('src/spring13/makecomparisonplot.R')
makeComparisonPlot(loadBed('data/s96rep1-high_peaks.bed'),loadBed('data/s96rep2-high_peaks.bed'),wiggleTable,4,6,'Comparison of raw read scores for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set1"))
mycor=cor(wiggleTable[,4],wiggleTable[,6])
text(50,220,substitute(paste(rho,"=",mycor),list(mycor=mycor)))


makeComparisonPlot(loadBed('data/s96rep1-high_peaks.bed'),loadBed('data/hs959rep1-new_peaks.bed'),wiggleTable,4,8,'Comparison of raw read scores for S96 vs HS959','S96', 'HS959', c("Overlap","S96 unique","HS959 unique"),brewer.pal(3,"Dark2"))

makeComparisonPlot(loadBed('data/s96rep1-high_peaks.bed'),loadBed('data/hs959rep1-new_peaks.bed'),backSubTable,4,8,'Comparison of raw read scores for S96 vs HS959','S96', 'HS959', c("Overlap","S96 unique","HS959 unique"),brewer.pal(3,"Spectral"))


mycor=cor(wiggleTable[,4],wiggleTable[,8])
text(40,250,substitute(paste(rho,"=",mycor),list(mycor=mycor)))
mypal<-brewer.pal(3,"Accent")







makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('s96rep2-high_peaks.bed'),normDiffTable,3,4,'Comparison of NormDiff scores for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set2"))
mycor=cor(normDiffTable[,3],normDiffTable[,4])
text(2,30,substitute(paste(rho,"=",mycor),list(mycor=mycor)))


makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('hs959rep1-new_peaks.bed'),normDiffTable,3,5,'Comparison of NormDiff scores for S96 vs HS959','S96','HS959',c("Overlap","S96","HS959"),brewer.pal(3,"BrBG"))
mycor=cor(normDiffTable[,4],normDiffTable[,8])
text(2,30,substitute(paste(rho,"=",mycor),list(mycor=mycor)))




makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('s96rep2-high_peaks.bed'),backSubTable,3,4,'Comparison of BackSub scores for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set2"))


makeComparisonPlotHelp(loadBed('s96overlap-high-peaks.bed'),loadBed('hs959rep1-new_peaks.bed'),backSubTable,3,5,'Comparison of BackSub scores for S96 vs HS959','S96','HS959',c("Overlap","S96","HS959"),brewer.pal(3,"Spectral"))



######################3
# MACS DIFF
######################3
makeComparisonPlotHelp(loadBed('s96overlap-high-peaks.bed'),loadBed('S96rep1-Diff_peaks.bed'),loadBed('S96rep2-Diff_peaks.bed'),wiggleTable,4,6,'Comparison of S96 replicates','S96rep1','S96rep2',c("Overlap","S96rep1","S96rep2"),brewer.pal(3,"Spectral"))


makeComparisonPlotHelp(loadBed('s96vshs959overlap.bed'),loadBed('S96vsHS959-Diff_peaks.bed'),loadBed('HS959vsS96-Diff_peaks.bed'),wiggleTable,4,8,'Comparison of S96 vs HS959','S96rep1','HS959rep1',c("Overlap","S96rep1","HS959rep1"),brewer.pal(3,"Spectral"))
#plot(retbin,main="Hexagonal binning")


chip<-read.BED('data/ELAND-BED/SEG1ChIP_repc.bed')
input<-read.BED('data/ELAND-BED/SEG1input.bed')
ret<-bin.data(chip,input,1000)
ret$input[ret$input>1000]<-1000
ret$chip[ret$chip>8000]<-8000
hexbin1<-hexbin(ret$input,ret$chip,xbins=250)
plot(hexbin1,colramp=rainbow)
