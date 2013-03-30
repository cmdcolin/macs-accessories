#Bedfile operations



loadBed<-function(path,header=FALSE) {
  bedfile<-read.table(path,sep='\t',header=header)
  
  names(bedfile)=c('chr','start','end','name','score')
  #bedfile<-transform(bedfile, start=as.numeric(start), end=as.numeric(end))
  bedfile
}


# revert back to old apply method of intersect rows
intersectBedLimma<-function(nc1,nc2) {
  chrsplit<-split(nc2,factor(nc2[,'chr']))
  x<-names(chrsplit)
  rteno<-sapply(strsplit(x,'\\.'),function(x)x[1])
  names(chrsplit)<-rteno
  sel<-apply(nc1, 1, function (row1) {
    rteno2<-strsplit(row1[['chr']],"\\.")[[1]][1]
    chrselect=chrsplit[[rteno2]]
    ret<-apply(chrselect, 1,function(row2){
      #X2 >= Y1 and Y2 >= X1
      r1<-as.numeric(row1['start'])
      r2<-as.numeric(row1['end'])
      s1<-as.numeric(row2['start'])
      s2<-as.numeric(row2['end'])
      (r2>s1)&&(s2>r1)
    })
    if(debug) {
      printf("%s %s %s %s --- %d,%d\n", row1['name'], row1['chr'], row1['start'], 
             row1['end'], nrow(chrselect),sum(ret))
    }

    sum(ret)>0
  })
  sel
}

uniqueBedLimma<-function(nc1,nc2) {
  chrsplit<-split(nc2,factor(nc2[,'chr']))
  x<-names(chrsplit)
  rteno<-sapply(strsplit(x,'\\.'),function(x)x[1])
  names(chrsplit)<-rteno
  sel<-apply(nc1, 1, function (row1) {
    rteno2<-strsplit(row1[['chr']],"\\.")[[1]][1]
    chrselect=chrsplit[[rteno2]]
    
    ret<-apply(chrselect, 1,function(row2){
      #  overlap AR < BL || BR < AL
      # Linear overlap AR < BL || BR < AL
      #X2 >= Y1 and Y2 >= X1
      #(y['start']<=x['end'])&&(x['start']<=y['end'])
      r1<-as.numeric(row1['start'])
      r2<-as.numeric(row1['end'])
      s1<-as.numeric(row2['start'])
      s2<-as.numeric(row2['end'])
      (r2<s1)||(s2<r1)
    })
    if(debug) {
      printf("%s %s %s %s --- %d,%d\n", row1['name'], row1['chr'], row1['start'], 
             row1['end'], nrow(chrselect),sum(ret))
    }
    
    sum(ret)==nrow(chrselect)
  })
  sel
}




# wrapper to return only interesting rows
intersectBed<-function(bed1,bed2) {
  bed1[intersectBedLimma(bed1,bed2),]
}


# wrapper to return only interesting rows
uniqueBed<-function(bed1,bed2) {
  bed1[uniqueBedLimma(bed1,bed2),]
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

