# library(ChIPpeakAnno)
# library(VennDiagram)
library(limma)
library(edgeR)
load('wiggleTable.RData')
# bed1=read.table('SEG1rep1u_peaks.bed')
# bed2=read.table('SEG1rep2u_peaks.bed')
# bed3=read.table('SEG1rep3u_peaks.bed')
# bed123=read.table('SEG1rep123_peaks.bed')
# 
# 
# bed1rd <- BED2RangedData(bed1,header=FALSE)
# bed2rd <- BED2RangedData(bed2,header=FALSE)
# bed3rd <- BED2RangedData(bed3,header=FALSE)
# bed123rd <- BED2RangedData(bed123,header=FALSE)



y=DGEList(wiggleTable[,c(-1,-2)],group=c(1,1,1,2,2,2),lib.size=c(2218566,1292603,2334845,1106212,1300837,2538686))
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
save(et,file="et.RData")





