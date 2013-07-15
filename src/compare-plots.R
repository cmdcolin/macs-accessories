################
# Main file for analyzezscore-joinmod.R 
#




#Comparison plots
source('src/spring13/makecomparisonplot.R')
makeComparisonPlot(loadBed('data/s96rep1-high_peaks.bed'),loadBed('data/s96rep2-high_peaks.bed'),wiggleTableScale,3,4,'Comparison of read counts for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set1"))
makeComparisonPlot(loadBed('data/s96rep1-high_peaks.bed'),loadBed('data/s96rep2-high_peaks.bed'),wiggleTableScale,3,4,'Comparison of  read counts for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(4,"Set1"),logscale=TRUE)
mycor=cor(wiggleTable[,4],wiggleTable[,6])
text(50,220,substitute(paste(rho,"=",mycor),list(mycor=mycor)))

makeComparisonPlot(loadBed('data/s96rep1-high_peaks.bed'),loadBed('data/hs959rep1-new_peaks.bed'),wiggleTableScale,3,6,'Comparison of read counts for S96 vs HS959','S96rep1', 'HS959rep1', c("Overlap","S96 unique","HS959 unique"),brewer.pal(3,"Dark2"))
makeComparisonPlot(loadBed('data/s96rep1-high_peaks.bed'),loadBed('data/hs959rep1-new_peaks.bed'),wiggleTableScale,3,6,'Comparison of read counts for S96 vs HS959','S96rep1', 'HS959rep1', c("Overlap","S96 unique","HS959 unique"),brewer.pal(4,"Dark2"),logscale=TRUE)




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
