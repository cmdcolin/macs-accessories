# Regex library
library(stringr)

#################
#! Constructor
WiggleClass<-function(name) {
  nc=list(
    name=name,
    scaling=1,
    variance=1,
    spacing=10,
    controlpath=paste(name,'/',name,'_MACS_wiggle/control/',sep=''),
    treatpath=paste(name,'/',name,'_MACS_wiggle/treat/',sep=''),
    controlname=paste(name,'_control_afterfiting_',sep=''),
    treatname=paste(name,'_treat_afterfiting_',sep=''),
    peaks=read.table(paste(name,'/',name,'_peaks.bed',sep='')),
    wiglist=list()
    )
  
  ########################
  # Read wiggle files from path into memory and assign filesnames
  nc$loadWiggles=function(e=environment()) {
    loadWiggle<-function(wigpath,env) {
      files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
      for (filename in files) {
        fn<-paste(wigpath,filename,sep='')
        if(debug)
          cat(filename, '\n')
        if(exists(filename,env))
          wig<-get(filename,env)
        else {
          wig<-read.table(fn, skip=2)
          assign(filename,wig,env)
        }
        nc$wiglist[[filename]]=wig
      }
    }
    
    loadWiggle(nc$treatpath,e)
    loadWiggle(nc$controlpath,e)
  }
  

  
  
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}



intersectBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    chr1=x[1]
    start1=as.integer(x[2])
    end1=as.integer(x[3])
    sublist=nc2[nc2$V1==chr1,]
    ret=apply(sublist,1,function(y){
      chr2=y[1]
      start2=as.integer(y[2])
      end2=as.integer(y[3])
      pn2=y[4]
      (start2<=end1)&&(start1<=end2)
    })
    sum(ret)>0
  })
  
  nc1[selectrows,]
}
uniqueBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    chr1=x[1]
    start1=as.integer(x[2])
    end1=as.integer(x[3])
    pn1=x[4]
    sublist=nc2[nc2$V1==chr1,]
    ret=apply(sublist,1,function(y){
      chr2=y[1]
      start2=as.integer(y[2])
      end2=as.integer(y[3])
      pn2=y[4]
      #AR < BL || BR < AL
      (start2<=end1)&&(start1<=end2)
    })
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

