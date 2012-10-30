# Sample average normdiff score from bedfile
#Make sure chromosome names match (i.e. see convert.R)


getMaxAvgZscoreAll2<-function(normfile,bedfile,ws=100) {  
  nd<-read.table(normfile,skip=1)
  bed<-read.table(bedfile)
  
  print(str(nd))
  print(str(bed))
  
  apply(bed,1,function(y) {
    chr=y[1]
    start=y[2]
    end=y[3]
    name=y[4]
    score=y[5]
    select_temp<-nd[y[1]==nd$V2,]
    test<-apply(select_temp,1,function(x){x[3]<=y[3]&&x[4]>=y[2]})
    #print(test)
    select<-select_temp[test,]
    mean(select$V5)
    #if(nrow(select)>0) {
    #  ret<-sapply(seq(1,nrow(select),by=ws),function(i) {
    #    start=i
    #    end=i+ws
    #    mean(select[start:end,4])
    #  })
    #  max(ret)
    #}
    #else {
    #  0
    #}
  })
}


f<-getwd()
setwd('data')
files=list.files(pattern="_normdiff.txt")
exp_names=str_replace_all(files,"_normdiff.txt","")
print(exp_names)
peaks="S96rep1_peaks.bed"
ret<-lapply(exp_names,function(name) {
  getMaxAvgZscoreAll2(sprintf("%s_normdiff.txt",name),peaks)
})


#Write output file
bed<-read.table(peaks)
x<-data.frame(name=bed$V4)
for(i in ret) {
  x<-cbind(x,i)
}
outfile<-x
names(outfile)<-c('YORF',exp_names)
names(outfile)
write.table(outfile,"test_output.txt",quote=FALSE,sep='\t',row.names=FALSE)
setwd(f)







################
# Old code


getMaxAvgZscoreAll<-function(wza,ws=100) {
  ret=data.frame(chr=character(0),start=numeric(0),end=numeric(0),normdiff=numeric(0))
  for(z in wza) {
    print
    chr=z[,1]
    tpos=z[,2]
    cpos=z[,3]
    scores=z[,4]
    if(debug)
      cat(length(scores),'\n')
    for(i in seq(1,length(scores)-1,by=ws)) {
      start=i
      end=i+ws
      mchr=chr[start]
      mscore=mean(as.numeric(scores[start:end]))
      row=data.frame(mchr,start,end,mscore)
      ret<-rbind(ret,row)
    }
    if(debug)
      cat(length(scores),' ',nrow(ret),'\n')
  }
  ret
}


getMaxAvgZscore<-function(wz,ws=10) {
  sapply(wz,function(zlist){
    
    zlist=unlist(zlist)
    len=length(zlist)
    reads=numeric(len)
    #if(debug)
    #  print(len)
    if(is.numeric(zlist)==FALSE) {
      cat('here1 ', length(zlist), typeof(zlist),'\n')
      #  print(unlist(zlist))
      #  zlist=unlist(zlist)
    }
    
    for(i in 1:len) {
      b=i
      e=i+ws
      reads[i]=mean(zlist[b:e],na.rm=TRUE)
    }
    max(reads)
  })
}