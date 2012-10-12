# Use regex library
library(stringr)
library(R.utils)


loadWiggle<-function(wigpath) {
  files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
  ret<-lapply(files,function(x){
    fn<-paste(wigpath,filename,sep='')
    if(debug)
      cat(filename, '\n')
    if(exists(filename))
      wig<-get(filename)
    else {
      wig<-read.table(fn, skip=2)
      assign(filename,wig)
    }
  })
  names(ret)<-files
  ret
}



setClass("WiggleClass",
         representation(
           name="character", 
           controlpath="character",
           treatpath="character",
           controlname="character",
           treatname="character",
           peaks="data.frame")
)
setMethod("initialize","WiggleClass",function(.Object,...,name) {
  printf("Loading %s\n", name)
  controlpath=sprintf("%s/%s_MACS_wiggle/control/",name,name)
  treatpath=sprintf("%s/%s_MACS_wiggle/treat/",name,name)
  controlname=sprintf("%s/%s_control_afterfiting_",name,name)
  treatname=sprintf("%s/%s_treat_afterfiting_",name,name)
  peaks=read.table(sprintf("%s/%s_peaks.bed",name,name))
  callNextMethod(.Object,...,name=name,
    controlpath=controlpath,treatpath=treatpath,
    controlname=controlname,treatname=treatname,
    peaks=peaks)
})

# Where are my ROMAN NUMERAL FUNCTIONS??
# Can't upload anything to RESEARCH COMPUTING??
# Evaluate MANORM pairwise on replicates
# Concatenate replicate data into one file
# Compare results
# Fix many roman numeral sequences

romanNum<-function(str) {
  if(str=="01") "I"
  else if(str=="02") "II"
  else if(str=="03") "III"
  else if(str=="04") "IV"
  else if(str=="05") "V"
  else if(str=="06") "VI"
  else if(str=="07") "VII"
  else if(str=="08") "VIII"
  else if(str=="09") "IX"
  else if(str=="10") "X"
  else if(str=="11") "XI"
  else if(str=="12") "XII"
  else if(str=="13") "XIII"
  else if(str=="14") "XIV"
  else if(str=="15") "XV"
  else if(str=="16") "XVI"
  else if(str=="mt") "M"
  else "Error"
}


convertFile<-function(filename) {
  con<-file(filename)
  open(con)
  lines<-readLines(con,-1)
  outlines<-sapply(lines[1:100], function(cline) {
    x<-str_match(cline,"(.*)(chr)([0-9a-z]{2})(.fsa)(.*)")
    if(sum(is.na(x))==0) {
      ret<-sprintf("%s%s%s%s",x[2],x[3],romanNum(x[4]),x[6])
      cat(ret)
      ret
    }
    else
      cline
  })
  writeLines(outlines,filename)
  close(con)
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

