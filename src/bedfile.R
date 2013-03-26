#Bed



loadBed<-function(path) {
  bedfile<-read.table(path,sep='\t',colClasses=c("character","numeric","numeric","character","numeric"))
  
  names(bedfile)=c('chr','start','end','name','score')
  #bedfile<-transform(bedfile, start=as.numeric(start), end=as.numeric(end))
  bedfile
}


#advice taken from http://stackoverflow.com/questions/5357003/get-and-process-entire-row-in-ddply-in-a-function
intersectBed<-function(nc1,nc2) {
  nc1$row<-1:nrow(nc1)
  nc2$row<-1:nrow(nc2)
  chrsplit<-split(nc2,factor(nc2$chr))
  sel<-daply(nc1, .(row), function (row1) {
    chrselect=chrsplit[[row1$chr]]
    ret<-daply(chrselect, .(row),function(row2){
      #  overlap AR < BL || BR < AL
      #X2 >= Y1 and Y2 >= X1
      intersect1D(row1$start,row1$end,row2$start,row2$end)
    })
    printf("%s %s %d %d --- %d,%d\n", row1$name, row1$chr, row1$start, row1$end, nrow(chrselect),sum(ret))

    sum(ret)>0
  })
  nc1[sel,]
}




uniqueBed<-function(nc1,nc2) {
  
  nc1$row<-1:nrow(nc1)
  nc2$row<-1:nrow(nc2)
  chrsplit<-split(nc2,factor(nc2$chr))
  sel<-daply(nc1, .(row), function (row1) {
    chrselect=chrsplit[[row1$chr]]
    ret<-daply(chrselect, .(row),function(row2){
      #  overlap AR < BL || BR < AL
      # Linear overlap AR < BL || BR < AL
      #X2 >= Y1 and Y2 >= X1
      #(y['start']<=x['end'])&&(x['start']<=y['end'])
      (row1$end<row2$start)||(row2$end<row1$start)
    })
    printf("%s %s %d %d --- %d,%d\n", row1$name, row1$chr, row1$start, row1$end, nrow(chrselect),sum(ret))
    
    sum(ret)==nrow(chrselect)
  })
  nc1[sel,]
  
  
  
  
#   selectrows=apply(nc1,1,function(x){
#     
#     sublist=nc2[nc2$chr==x['chr'],]
#     ret=apply(sublist,1,function(y){
#       (y['start']<=x['end'])&&(x['start']<=y['end'])
#     })
#     
#     ##! Get unique peaks with no overlap
#     sum(ret)==0
#   })
#   nc1[selectrows,]
}



intersect1D<-function(a1,a2,b1,b2) {
  (a2>b1)&&(b2>a1)
}

###########


######Sketches
#
# Get peak index from bed
#  (uneeded) match(bed1,bed2)
###########
#getPeakIndex<-function(bedfile)
#{
#  index=array()
#  for(i in 1:length(bedfile$V1)) {
#    index[i]=as.integer(unlist(strsplit(as.character(bedfile$V4[i]),'_'))[3])
#  }
#  index
#}
#########
#
#nc$overlapBed(b) {
#  s=paste('intersectBed -a ',nc$name,'/',nc$name,'_peaks.bed')
#  s=paste(s, '-b ', b$name, '/', b$name, '_peaks.bed -wa > ')
#  s=paste(nc$name, '/', nc$name,'_overlap.bed')
#  system(s)
#}
#nc$uniqueBed(b) {
#  s=paste('subtractBed -a ',nc$name,'/',nc$name,'_peaks.bed')
#  s=paste(s, '-b ', nc$name, '/', nc$name, '_overlap.bed > ')
#  s=paste(a$name, '/', a$name,'_unique.bed')
#  system(s)
#}
#plotnew(x1,x2,b1,b2,b3,name1,name2,func) {
#  id1=match(b2$V4,b1$V4)
#  id2=match(b3$V4,b1$V4)
#  r1=hs959$func(b1)
#  r2=s96$func(b1)
#  
#}

