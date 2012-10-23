

####
# Use all chromosomes for scaling factor

estimateScalingFactor <- function(wig) {
  sfactor<-function(x) {
    if(debug) {
      print(names(x))
    }
    treat=x[1]
    control=x[2]
    corr=findInterval(treat[,1],control[,1])
    control[corr,2]/treat[corr,2]
  }
  treat=wig@wiggles[["treat"]]
  control=wig@wiggles[["control"]]
  ratio_data=apply(cbind(treat,control),1,sfactor)
  median(unlist(ratio_data),na.rm=TRUE)
}

# Sqrt(Aw+Bw/c), w=all
estimateVarianceAll<-function(wig,scaling) {
  getSignal<-function(file) {
    wig$wiglist[[file]]$V2
  }
  chip_signal=sapply(wig@wiggles[["treat"]],getSignal)
  control_signal=sapply(wig@wiggles[["control"]],getSignal)
  #Average signal
  average_chip=mean(chip_signal)
  average_control=mean(control_signal)
  sqrt(average_chip+average_control/scaling^2)
}


####
estimateVarianceWindow<-function(xpos, treat,control,ws,corr1,corr2,scaling) {
  # find start end positions
  # check for inconsistencies in data
  #
  # checkspacing is correct b window
  #select signals
  #cat(xpos,'\t(', treat$V1[beginning],',',treat$V1[ending],')',nc$matches, '-', nc$unmatches,'\t',
  #    nc$matches1,' ', nc$unmatches1,'(',beginning, ' ', ending, ') (', nc$shifter, ' ', nc$shifter2, '\n')
  
  b1=xpos[1]-ws
  e1=xpos[1]+ws
  b2=xpos[2]-ws
  e2=xpos[2]+ws
  b1=if(b1<1) 1 else b1
  b2=if(b2<1) 1 else b2
  e1=if(e1>dim(treat)[1]) dim(treat)[1] else e1
  e2=if(e2>dim(control)[1]) dim(control)[1] else e2
  
  chip_signal=treat$V2[b1:e1]
  control_signal=control$V2[b2:e2]
  average_chip=mean(chip_signal)
  average_control=mean(control_signal)
  sqrt(average_chip+average_control/scaling^2)
}


Zxi<-function(x, treat,control,window,corr1,corr2,scaling,variance) {
  pos1=x[1]
  pos2=x[2]
  ma=sapply(window,function(w){
    estimateVarianceWindow(x,treat,control,w,corr1,corr2,scaling)
  })
  
  ## Normdiff
  A=treat[pos1,2]
  B=control[pos2,2]
  c=scaling
  sigma=max(ma,variance)
  (A-B/c)/sigma
}



# Calculate Z scores over all wiggle files
Z<-function(wig, bedfile, scaling, variance, window=c(1,10)) {
  # Get max average reads over window size
  getZscore<-function(x,window){
    chr=x[1];
    start=as.integer(x[2])
    end=as.integer(x[3])
    tf=paste(wig$treatname,chr,'.wig.gz',sep='')
    cf=paste(wig$controlname,chr,'.wig.gz',sep='')
    treat=wig$wiglist[[tf]]
    control=wig$wiglist[[cf]]
    mw=max(window)
    corr1=findInterval(seq(start-mw,end+mw,by=nc$spacing),treat$V1)
    corr2=findInterval(seq(start-mw,end+mw,by=nc$spacing),control$V1)
    if(debug==TRUE)
      cat(chr,'-\t(',start,',',end, ')\n')
    app=cbind(corr1,corr2)
    apply(app,1,Zxi,treat,control,window,corr1,corr2)
  }
  
  apply(bedfile,1,getZscore,window)
}





# Calculate Z scores over all wiggle files
Zall<-function(wig, window=c(1,10)) {
  files1=list.files(path=nc$treatpath,pattern="*.fsa.wig.gz")
  files2=list.files(path=nc$controlpath,pattern="*.fsa.wig.gz")
  ret=apply(cbind(files1,files2),1,function(f){
    treat<-wig$wiglist[[f[1]]]
    control<-wig$wiglist[[f[2]]]
    if(debug==TRUE)
      cat(f[1],'\t',f[2],'\n')
    
    corr1=1:length(treat$V1)
    corr2=findInterval(treat$V1,control$V1,all.inside=TRUE)
    app=cbind(corr1,corr2)
    
    listret=apply(app, 1, function(x){nc$Zxi(x,treat,control,window,corr1,corr2)})
    
    # first column chr.fsa
    a=character(length(corr1))
    chr=str_extract(f[1],"chr[0-9a-z]{2}.fsa")
    for(i in 1:length(corr1)) {
      a[i]=chr
    }
    cbind(a,treat$V1[corr1],control$V1[corr2],listret)
  })
  
  # get chr w/ regex
  #todo get chr from filename
  #for(i in 1:length(ret)) {
  #  name=files1[i]
  #  chr=str_extract(name,"chr[0-9a-z]{2}.fsa")
  #  ret[i][['chr']]=chr
  #}
  #ret
}





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