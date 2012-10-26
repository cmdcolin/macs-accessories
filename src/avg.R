

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