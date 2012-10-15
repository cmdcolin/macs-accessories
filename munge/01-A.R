# Example preprocessing script.

source('src/wiggle.R')

#f<-getwd()
#setwd('data')
#dataset=list()
#for(filename in list.files()) {
#  dataset[filename]=new("WiggleClass",name=filename)
#}
#setwd(f)

#--kkkkkkk,,,,kkkj.;/;/;;''''''''''''''''''''



# Convert wiggle file chr names
# 
# for(file in c(list.files(recursive=TRUE,pattern="*.wig$"),
#               list.files(recursive=TRUE,pattern="*.bed$"))) 
# {
#   printf("Converting %s\n",file)
#   # Uses sacCer3 chromosome
#   if(saccer==TRUE)
#     convertFileSacCer(file)
#   # Uses S288c chromosome
#   else if(s288c==TRUE)
#     convertFileS288C(file)
# }


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
  print(x[1])
  first=l[[x[1]]]
  second=l[[x[2]]]
  inter=intersectBedLimma(first,second)
  nrow(inter)
})
# Make a venn diagram from all intersect overlap!!


ret<-intersectBed(x1,x2)
ret2<-uniqueBed(x1,x2)
ret3<-uniqueBed(x2,x1)

printf("Found %d overlapping peaks\n", nrow(ret))
printf("Found %d unique peaks in x1\n", nrow(ret2))
printf("Also found %d unique peaks in x2\n", nrow(ret3))
