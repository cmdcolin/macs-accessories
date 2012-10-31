# Test
spacing=10


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


estimateVarianceAllMod<-function(wig1,wig2,s1,s2) {
  files1=ls(wig1$treat)
  files2=ls(wig1$control)
  files3=ls(wig2$treat)
  files4=ls(wig2$control)
  getSignal<-function(file,curr) {
    curr[[file]][,2]
  }
  sig1=unlist(lapply(files1,getSignal,wig1$treat))
  sig2=unlist(lapply(files2,getSignal,wig1$control))
  sig3=unlist(lapply(files3,getSignal,wig2$treat))
  sig4=unlist(lapply(files4,getSignal,wig2$control))
  #Average signal
  average_chip=mean(sig1)
  average_control=mean(sig2)
  average_chip2=mean(sig3)
  average_control2=mean(sig4)
  sqrt(average_chip+average_control/s1^2+average_chip2+average_control2/s2^2)
}


estimateScalingFactor <- function(wig) {
  sfactor<-function(x) {
    if(debug) {
      #print(str(x))
    }
    t=wig$treat[[x[1]]]
    c=wig$control[[x[2]]]
    pos=findInterval(t[,1],c[,1])
    c[pos,2]/t[pos,2]
  }
  
  ratio_data=apply(cbind(ls(wig$treat),ls(wig$control)),1,sfactor)
  median(unlist(ratio_data),na.rm=TRUE)
}












plotAddSubNormDiffV2<-function(w1,w2,peakfile1,peakfile2) {
  peak1=loadBed(peakfile1)
  peak2=loadBed(peakfile2)
  shared=intersectBed(peak1,peak2)
  unique=uniqueBed(peak1,peak2)
  id1=match(shared$name,peak1$name) 
  id2=match(unique$name,peak1$name)
  s1=estimateScalingFactor(w1)
  s2=estimateScalingFactor(w2)
  vall1=estimateVarianceAllMod(w1,w2,s1,s2)
  ret1=Zmod(peak1, w1,w2,Zaddxi,vall1,s1,s2)
  ret2=Zmod(peak1, w1,w2,Zsubxi,vall1,s1,s2)

  ms1<-lapply(ret1,function(l) { 
    max(sapply(seq(1,length(l)-1,by=10),function(i) {
      start=i
      end=i+10
      mean(l[start:end],na.rm=TRUE)
    }))
    })
  
  ms2<-lapply(ret2,function(l) { 
    max(sapply(seq(1,length(l)-1,by=10),function(i) {
      start=i
      end=i+10
      mean(l[start:end],na.rm=TRUE)
    }))
  })
  r1=sapply(ret1,mean)
  r2=sapply(ret2,mean)
  #plot(r1,r2,pch='*',xlab='Additive NormDiff',ylab='Subtractive NormDiff')
  #points(r1[id1],r2[id1],col='blue')
  #points(r1[id2],r2[id2],col='red')
  
  
  plot(ms1,ms2,pch='*',xlab='Additive NormDiff',ylab='Subtractive NormDiff')
  points(ms1[id1],ms2[id1],col='blue')
  points(ms1[id2],ms2[id2],col='orange')
}