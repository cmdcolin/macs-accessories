#Bed



loadBed<-function(path) {
  bedfile<-read.table(path,sep='\t',colClasses=c("character","numeric","numeric","character","numeric"))
  
  names(bedfile)=c('chr','start','end','name','score')
  
  bedfile
}



intersectBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    
    sublist=nc2[nc2$chromosome==x['chromosome'],]
    ret=apply(sublist,1,function(y){
      #  overlap AR < BL || BR < AL
      (y['start']<=x['end'])&&(x['start']<=y['end'])
    })
    
    ##! Get overlap peaks where intersect>0
    sum(ret)>0
  })
  nc1[selectrows,]
}



intersectBedLimma<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    
    sublist=nc2[nc2$chromosome==x['chromosome'],]
    ret=apply(sublist,1,function(y){
      #  overlap AR < BL || BR < AL
      (y['start']<=x['end'])&&(x['start']<=y['end'])
    })
    
    ##! Get overlap peaks where intersect>0
    sum(ret)>0
  })
  selectrows
}


uniqueBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    
    sublist=nc2[nc2$chromosome==x['chromosome'],]
    ret=apply(sublist,1,function(y){
      # Linear overlap AR < BL || BR < AL
      (y['start']<=x['end'])&&(x['start']<=y['end'])
    })
    
    ##! Get unique peaks with no overlap
    sum(ret)==0
  })
  nc1[selectrows,]
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

