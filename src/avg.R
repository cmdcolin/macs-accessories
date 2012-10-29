getMaxAvgZscoreAll2<-function(normdiff,bedfile,ws=100) {  
  wza<-read.table(normdiff)
  bed<-read.table(bedfile)
  
  
  apply(bed,1,function(row) {
    chr=row[1]
    start=row[2]
    end=row[3]
    name=row[4]
    score=row[5]
    select<-wza[wza[,1]==chr && wza[,2]>=start && wza[,3]<=end,]
    
    
    ret<-lapply(seq(1,nrow(select),by=ws),function(i) {
      start=i
      end=i+ws
      mchr=chr[start]
      mscore=mean(wza[start:end,4])
    })
    max(ret)
  })
}


f<-getwd()
setwd('data')
files=list.files(pattern="_normdiff.txt")
names=str_replace_all(files,"_normdiff.txt","")
print(names)
ret<-lapply(names,function(name) {
  getMaxAvgZscoreAll2(sprintf("%s_normdiff.txt",name),sprintf("%s_peaks.bed",name))
})
setwd(f)


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