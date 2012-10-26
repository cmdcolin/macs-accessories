
library(gtools)
x1<-loadBed('data/HS959rep1_peaks.bed')
x2<-loadBed('data/HS959rep2_peaks.bed')
x3<-loadBed('data/S96rep1_peaks.bed')
x4<-loadBed('data/S96rep2_peaks.bed')
l=list()
l[[1]]=x1
l[[2]]=x2
l[[3]]=x3
l[[4]]=x4
ind<-permutations(4,2,1:4)
ret<-apply(ind,1,function(x) {
  printf("Combining %s\n",str(x))
  first=l[[x[1]]]
  second=l[[x[2]]]
  intersectBedLimma(first,second)
})
# Make a venn diagram from all intersect overlap!!
