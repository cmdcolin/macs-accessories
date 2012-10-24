spacing=10

####
# Use all chromosomes for scaling factor

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

# Sqrt(Aw+Bw/c), w=all
estimateVarianceAll<-function(wig,scaling) {
  getSignal<-function(file,curr) {
    curr[[file]][,2]
  }
  chip_signal=unlist(lapply(ls(wig$treat),getSignal,wig$treat))
  control_signal=unlist(lapply(ls(wig$control),getSignal,wig$control))
  
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
    treat=wig$treat[[chr]]
    control=wig$control[[chr]]
    mw=max(window)
    corr1=findInterval(seq(start-mw,end+mw,by=spacing),treat$V1)
    corr2=findInterval(seq(start-mw,end+mw,by=spacing),control$V1)
    if(debug==TRUE)
      cat(chr,'-\t(',start,',',end, ')\n')
    app=cbind(corr1,corr2)
    apply(app,1,Zxi,treat,control,window,corr1,corr2,scaling,variance)
  }
  
  apply(bedfile,1,getZscore,window)
}





# Calculate Z scores over all wiggle files
Zall<-function(wig, scaling, variance, window=c(1,10)) {
  files1=ls(wig$treat)
  files2=ls(wig$control)
  ret=apply(cbind(files1,files2),1,function(f){
    t<-wig$treat[[f[1]]]
    c<-wig$control[[f[2]]]
    if(debug==TRUE)
      cat(f[1],'\t',f[2],'\n')
    
    corr1=1:length(t$V1)
    corr2=findInterval(t$V1,c$V1,all.inside=TRUE)
    app=cbind(corr1,corr2)
    
    listret=apply(app, 1, function(x){Zxi(x,t,c,window,corr1,corr2,scaling,variance)})
    
    # first column chr.fsa
    chr=f[1]
    cbind(rep(chr,length(corr1)),t$V1[corr1],c$V1[corr2],listret)

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



