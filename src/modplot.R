# Test

# Calculate Z scores over all wiggle files
Zmod<-function(bedfile, wig1, wig2, func, vall, s1,s2,window=c(1,10)) {
  # Get max average reads over window size
  getZscore<-function(x,wig1,wig2){
    chr=x[1];
    start=as.integer(x[2])
    end=as.integer(x[3])
    mw=max(window)
    
    # Get wig1 corr
    treat1=wig1$treat[[chr]]
    control1=wig1$control[[chr]]
    corr1=findInterval(seq(start-mw,end+mw,by=spacing),treat1$V1)
    corr2=findInterval(seq(start-mw,end+mw,by=spacing),control1$V1)
    
    # Get wig2 corr
    treat2=wig2$treat[[chr]]
    control2=wig2$control[[chr]]
    corr3=findInterval(seq(start-mw,end+mw,by=spacing),treat2$V1)
    corr4=findInterval(seq(start-mw,end+mw,by=spacing),control2$V1)

    if(debug==TRUE)
      cat(chr,'-\t(',start,',',end, ')\n')
    app=cbind(corr1,corr2,corr3,corr4)
    apply(app,1,func,treat1,treat2,control1,control2,window,corr1,corr2,corr3,corr4,vall,s1,s2)
  }
  
  apply(bedfile,1,getZscore,wig1,wig2)
}


Zaddxi<-function(x,t1,t2,c1,c2,window,corr1,corr2,corr3,corr4,vall,s1,s2) {
  pos1=x[1]
  pos2=x[2]
  pos3=x[3]
  pos4=x[4]
  ma=sapply(window,function(w){
    estimateVarianceWindowMod(x,t1,t2,c1,c2,w,corr1,corr2,corr3,corr4,s1,s2)
  })
  
  ## Normdiff
  A1=t1[pos1,2]
  B1=c1[pos2,2]
  A2=t2[pos3,2]
  B2=c2[pos4,2]
  sigma=max(ma,vall)
  ret=((A1-B2/s1)+(A2-B2/s2))/sigma
  ret
}

Zsubxi<-function(x,t1,t2,c1,c2,window,corr1,corr2,corr3,corr4,vall,s1,s2) {
  pos1=x[1]
  pos2=x[2]
  pos3=x[3]
  pos4=x[4]
  ma=sapply(window,function(w){
    estimateVarianceWindowMod(x,t1,t2,c1,c2,w,corr1,corr2,corr3,corr4,s1,s2)
  })
  
  ## Normdiff
  A1=t1[pos1,2]
  B1=c1[pos2,2]
  A2=t2[pos3,2]
  B2=c2[pos4,2]
  c1=s1
  c2=s2
  sigma=max(ma,vall)
  ret=((A1-B2/c1)-(A2-B2/c2))/sigma
  ret
}



####
estimateVarianceWindowMod<-function(xpos, t1,t2,c1,c2,ws,corr1,corr2,corr3,corr4,s1,s2) {
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
  b3=xpos[3]-ws
  e3=xpos[3]+ws
  b4=xpos[4]-ws
  e4=xpos[4]+ws
  b1=if(b1<1) 1 else b1
  b2=if(b2<1) 1 else b2
  b3=if(b3<1) 1 else b3
  b4=if(b4<1) 1 else b4
  e1=if(e1>nrow(t1)) nrow(t1) else e1
  e2=if(e2>nrow(c1)) nrow(c1) else e2
  e3=if(e3>nrow(t2)) nrow(t2) else e3
  e4=if(e4>nrow(c2)) nrow(c2) else e4
  
  chip_signal=t1$V2[b1:e1]
  control_signal=c1$V2[b2:e2]
  chip_signal2=t2$V2[b3:e3]
  control_signal2=c2$V2[b4:e4]
  average_chip=mean(chip_signal)
  average_control=mean(control_signal)
  average_chip2=mean(chip_signal2)
  average_control2=mean(control_signal2)
  sqrt(average_chip+average_control/s1^2+average_chip2+average_control2/s2^2)
}


estimateVarianceAllMod<-function(wig1,wig2) {
  files1=list.files(wig1$treatpath,pattern="*.fsa.wig.gz")
  files2=list.files(wig1$controlpath,pattern="*.fsa.wig.gz")
  files3=list.files(wig2$treatpath,pattern="*.fsa.wig.gz")
  files4=list.files(wig2$controlpath,pattern="*.fsa.wig.gz")
  getSignal<-function(file,wig){wig$wiglist[[file]]$V2}
  sig1=unlist(lapply(files1,getSignal,wig1))
  sig2=unlist(lapply(files2,getSignal,wig1))
  sig3=unlist(lapply(files3,getSignal,wig2))
  sig4=unlist(lapply(files4,getSignal,wig2))
  #Average signal
  average_chip=mean(sig1)
  average_control=mean(sig2)
  average_chip2=mean(sig3)
  average_control2=mean(sig4)
  sqrt(average_chip+average_control/wig1$scaling^2+average_chip2+average_control2/wig2$scaling^2)
}
