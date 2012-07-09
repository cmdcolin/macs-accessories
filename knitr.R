
plotTotalReads<-function(t, w1, w2, c1, c2) {
  r1=w1$getTotalReads(w1$peaks) 
  r2=w2$getTotalReads(w1$peaks)
  plot(r1,r2,xlab=paste(w1$name,'reads'),ylab=paste(w2$name, 'reads'),pch='*') 
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  points(r1[id1],r2[id1],pch=1,col=c1) 
  points(r1[id2],r2[id2],pch=1,col=c2) 
  title(t)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
}




plotMaxAvgReads<-function(t, w1, w2,c1,c2) {
  ###########
  # Use Max avg reads over windows
  rma1=w1$getMaxAvgReads(w1$peaks)
  rma2=w2$getMaxAvgReads(w1$peaks)
  ###################### 
  plot(rma1,rma2,xlab=paste(w1$name,'peak reads'),ylab=paste(w2$name,'peak reads'),pch='*')
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  points(rma1[id1],rma2[id1],pch=1,col=c1)
  points(rma1[id2],rma2[id2],pch=1,col=c2)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}



plotAvgZscore<-function(t, w1, w2, wz1, wz2, c1,c2) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  Z1<-sapply(wz1,mean)
  Z2<-sapply(wz2,mean)
  plot(Z1,Z2,xlab=paste('Max Avg', w1$name,'Zscore'),ylab=paste('Max Avg',w2$name,'Zscore'),pch='*')
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  points(Z1[id1],Z2[id1],pch=1,col=c1)
  points(Z1[id2],Z2[id2],pch=1,col=c2)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
  ret=list()
  ret[['1shared']]=Z1[id1]
  ret[['1unique']]=Z1[id2]
  ret[['2shared']]=Z2[id1]
  ret[['2unique']]=Z2[id2]
  ret
}

plotMaxAvgZscore<-function(t, w1, w2, wz1,wz2, c1,c2) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  maxw1<-w1$getMaxAvgZscore(wz1)
  maxw2<-w2$getMaxAvgZscore(wz2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  points(maxw1[id1],maxw2[id1],col=c1)
  points(maxw1[id2],maxw2[id2],col=c2)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
  ret=list()
  ret[['1shared']]=maxw1[id1]
  ret[['1unique']]=maxw1[id2]
  ret[['2shared']]=maxw2[id1]
  ret[['2unique']]=maxw2[id2]
  ret
}

plotSortedZscoreLines<-function(t, w1, w2, r, c1,c2,c3,c4) {
  wz1sort=sort(c(r[['2shared']],r[['2unique']]))
  wz2sort=sort(c(r[['1shared']],r[['1unique']]))
  xs1=1:length(wz1sort)
  xs11=match(r[['2unique']],wz1sort)
  xs12=match(r[['2shared']],wz1sort)
  xs2=1:length(wz2sort)
  xs21=match(r[['1unique']],wz2sort)
  xs22=match(r[['1shared']],wz2sort)
  plot(xs1,wz1sort,pch='.',xlab='Rank',ylab='Avg NormDiff')
  points(xs12,wz1sort[xs12],col=c3,pch='.')
  points(xs11,wz1sort[xs11],col=c4,pch='.')
  points(xs21,wz2sort[xs21],col=c1,pch='.')
  points(xs22,wz2sort[xs22],col=c2,pch='.')
  for(i in 1:length(xs12))
    lines(c(xs12[i],xs22[i]),c(wz1sort[xs12][i],wz2sort[xs22][i]),col=c1)
  for(i in 1:length(xs11))
    lines(c(xs11[i],xs21[i]),c(wz1sort[xs11][i],wz2sort[xs21][i]),col=c2)
  #polygon(c(xs1,rev(xs)),c(wz1sort,rev(wz2sort)),col='lightyellow',border=FALSE)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c3, c4))
  title(t)
}

plotSortedMaxAvgZscoreX<-function(t, w1, w2, r, c1,c2) {
  wz1sort=sort(r[['2shared']])
  wz2sort=sort(r[['2unique']])
  xs1=1:length(wz1sort)
  xs2=1:length(wz2sort)
  plot(xs1,wz1sort,pch='.',xlab='Rank',ylab='Avg NormDiff')
  points(xs1,wz1sort,col=c1,pch='.')
  points(xs2,wz2sort,col=c2,pch='.')
  #polygon(c(xs1,rev(xs)),c(wz1sort,rev(wz2sort)),col='lightyellow',border=FALSE)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}




plotZall2<-function(cols,w1,w2){
  select=cols[order(cols[,1]),]
  plot(x=1:nrow(select),y=select$mscore,pch='.',xlab='Sort rank', ylab='NormDiff score')
  for(i in 1:nrow(select)) {
    mycol=select$color[i]
    mypch='.'
    if(mycol!="yellow")
      mypch='*'
    points(x=i,y=select$mscore[i],col=mycol,pch=mypch)
  }
  legend("topleft",
         legend=c(paste('peak in',w2$name),paste('peak in',w1$name),'peak in both'),
         title='window overlaps',
         fill=c('blue','green','red'))
  title(paste('Sorted normdiff scores for entire',w1$name,'genome w=100'))
}
#Defunct
plotZall3<-function(cols,w1,w2){
  sel=cols[order(cols[,1]),]
  s1=sel[sel$color=="blue",]
  s2=sel[sel$color=="red",]
  s3=sel[sel$color=="green",]
  plot(1:nrow(sel),sel$mscore,pch='.')
  points(1:nrow(s3),s3$mscore,col="green")
  points(1:nrow(s1),s1$mscore,col="blue")
  points(1:nrow(s2),s2$mscore,col="red")
  legend("topleft",
         legend=c(paste('peak in',w2$name),paste('peak in',w1$name),'peak in both'),
         fill=c('blue','green','red'))
  title(paste('Sorted normdiff scores for entire',w1$name,'genome w=100'))
}


plotZall<-function(ts,w1,w2,cols=NULL) {
  if(is.null(cols)) {
    cols=data.frame(score=numeric(0),color=character(0))
    
    for(i in 1:nrow(ts)) {
      z=ts[i,]
      if(debug)
        print(z)
      #l=as.table(cbind(z[1],z[2],z[3]))
      a=intersectBed(z,w1$peaks)
      b=intersectBed(z,w2$peaks)
      if(nrow(b)!=0&&nrow(a)!=0)colselect='red'
      else if(nrow(b)!=0) colselect='blue'
      else if(nrow(a)!=0)colselect='green'
      else colselect='yellow'
      cols<-rbind(cols,data.frame(score=z[4],color=colselect))
    }
  }
  cols
  
  
##  plot(zaw)
}




plotMaxAvgZscore<-function(t, w1, w2, z1,z2, c1,c2) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  points(maxw1[id1],maxw2[id1],col=c1)
  points(maxw1[id2],maxw2[id2],col=c2)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}




plotZscoreCutoff<-function(t, w1, w2, z1,z2, cutoff) {

  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  
  # get subsets
  scales=1-pnorm(maxw2)  
  set1=(scales[id2]<cutoff)
  set2=(scales[id2]>=cutoff)
  # plot colors
  points(maxw1[id1],maxw2[id1],col='blue')
  points(maxw1[id2][set2],maxw2[id2][set2],col='red')
  points(maxw1[id2][set1],maxw2[id2][set1],col='yellow',pch=20)
  legend('bottomright', legend=c('shared', 'unique ', 'new'), fill=c('blue', 'red','yellow'))
  
  title(t)
  ret=scales<cutoff
  ret[id1]=FALSE
  ret
}




plotZscoreCutoffShared<-function(t, w1, w2, z1,z2, cutoff) {
  
  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  
  # get subsets
  scales=1-pnorm(maxw1)  
  set1=(scales[id1]<cutoff)
  set2=(scales[id1]>=cutoff)
  # plot colors
  points(maxw1[id2],maxw2[id2],col='blue')
  points(maxw1[id1][set2],maxw2[id1][set2],col='red')
  points(maxw1[id1][set1],maxw2[id1][set1],col='yellow',pch=20)
  legend('bottomright', legend=c('shared', 'unique ', 'new'), fill=c('blue', 'red','yellow'))
  
  title(t)
  ret=scales<cutoff
  #ret[id1]=FALSE
  ret
}




plotZscoreColor<-function(t, w1, w2, z1,z2) {

  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  
  scales=1-pnorm(maxw2)
  hsvscale=-log10(1-pnorm(maxw2))
  maxv=4/3*max(hsvscale)
  hsvscale=hsvscale/maxv
  hsvscale[hsvscale>1]=0.75
  for(i in 1:length(maxw2)) {
    points(maxw1[i],maxw2[i],col=hsv(hsvscale[i]))
  }
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}





plotZscoreMACS<-function(t, w1, w2, z1,z2, cutoff) {

  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  
  
  scales=1-pnorm(maxw1)
  hsvscale=-log10(1-pnorm(maxw1))
  #maxv=4/3*max(hsvscale)
  hsvscale=hsvscale/maxv
  hsvscale[hsvscale>1]=0.75
  
  # Get MACS pvalues
  tab=read.table(paste(wig1$name,'/',wig1$name,'_peaks.xls',sep=''),sep='\t',skip=2)
  pvals=as.numeric(as.character(tab$V7[2:nrow(tab)]))

  # Use MACS pavvls
  maxv=4/3*max(pvals)
  hsvscale=pvals/maxv
  hsvscale[hsvscale>1]=3/4
  for(i in 1:length(maxw1)) {
      points(maxw1[i],maxw2[i],col=hsv(hsvscale[i]))
  }
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}



## GGPlot2
plotGGPlot2<-function(t, w1, w2, z1,z2, cutoff) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  p=qplot(maxw1,maxw2,xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  p<-p+ scale_color_hue()
}








plotMAPlot2<-function(t, w1, w2, z1,z2, z3,z4,cutoff) {

  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  maxw3<-w1$getMaxAvgZscore(z3)
  maxw4<-w2$getMaxAvgZscore(z4)
  
  # shared /unique
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  shared=intersectBed(w2$peaks,w1$peaks)
  unique=uniqueBed(w2$peaks,w1$peaks)
  id3=match(shared$V4,w2$peaks$V4) 
  id4=match(unique$V4,w2$peaks$V4)
  
  # MA plot scores
  M1=log2(maxw1)-log2(maxw2)
  A1=1/2*(log2(maxw1)+log2(maxw2))
  M2=log2(maxw3)-log2(maxw4)
  A2=1/2*(log2(maxw3)+log2(maxw4))
  
  plot(A1,M1,pch='*',ylim=c(-6,6))
  points(A2,M2,pch='*')
  
  points(A1[id1],M1[id1],col='blue')
  points(A1[id2],M1[id2],col='green')
  points(A2[id3],M2[id3],col='blue')
  points(A2[id4],M2[id4],col='red')
  
  title(t)
}



plotMAPlot1<-function(t, w1, w2, z1,z2, z3,z4,cutoff) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  maxw3<-w1$getMaxAvgZscore(z3)
  maxw4<-w2$getMaxAvgZscore(z4)
  
  
  M1=log2(maxw1)-log2(maxw2)
  A1=1/2*(log2(maxw1)+log2(maxw2))
  
  M2=log2(maxw3)-log2(maxw4)
  A2=1/2*(log2(maxw3)+log2(maxw4))
  
  plot(A1,M1,pch='*',ylim=c(-6,6))
  points(A2,M2,pch='*')
  
  scales1=1-pnorm(maxw2)
  scales2=1-pnorm(maxw3)
  s1=-log10(scales1)
  s2=-log10(scales2)
  smax=4/3*max(s1,s2)
  for(i in 1:length(maxw1)) {
  
    points(A1[i],M1[i],col=hsv(s1[i]/smax,alpha='0.6'))
  }  
  for(i in 1:length(maxw3)) {
    
    points(A2[i],M2[i],col=hsv(s2[i]/smax,alpha='0.6'))
  }
  #scales1=1-pnorm(maxw1)
  #scales2=1-pnorm(maxw2)
  #s1=-log10(scales1)
  #s2=-log10(scales2)
  #plot(s1,s2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  
  #points(s1[scales2<cutoff],s2[scales2<cutoff],col='green')
  #points(s1[scales2>=cutoff],s2[scales2>=cutoff],col='red')
  
  title(t)
}


plotMaxAvgZscoreW<-function(t, w1, w2, z1,z2, cutoff) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  
  
  
  scales1=1-pnorm(maxw1)
  scales2=1-pnorm(maxw2)
  s1=-log10(scales1)
  s2=-log10(scales2)
  plot(s1,s2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  
  points(s1[scales2<cutoff],s2[scales2<cutoff],col='green')
  points(s1[scales2>=cutoff],s2[scales2>=cutoff],col='red')
  
  title()
}




plotMaxAvgZscoreZ<-function(t, w1, w2, z1,z2, cutoff) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  maxw1<-w1$getMaxAvgZscore(z1)
  maxw2<-w2$getMaxAvgZscore(z2)
  shared=intersectBed(w1$peaks,w2$peaks)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id1=match(shared$V4,w1$peaks$V4) 
  id2=match(unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  
  scales=1-pnorm(maxw2)
  hsvscale=-log10(1-pnorm(maxw2))
  maxv=3/4*max(hsvscale)
  hsvscale=hsvscale/maxv
  hsvscale[hsvscale>1]=0.75
  id3=match(scales[scales[id2]<cutoff],scales)
  for(i in 1:length(id1)) {
    points(maxw1[id1][i],maxw2[id1][i],col='green')
  }
  for(i in 1:length(id2)) {
    if(scales[id2][i]<cutoff) {
      points(maxw1[id2][i],maxw2[id2][i],
             col=hsv(hsvscale[id2][i]),
             pch=20)
    }
    else points(maxw1[id2][i],maxw2[id2][i],col='red')
  }
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}





getZcorrelation<-function(wz1,wz2) {
  apply(cbind(wz1,wz2),1,function(z){
    z1=z[1]
    z2=z[2]
    cor(unlist(z1),unlist(z2))
  })
}



getZcovariance<-function(wz1,wz2) {
  apply(cbind(wz1,wz2),1,function(z){
    z1=z[1]
    z2=z[2]
    cov(sort(unlist(z1)),sort(unlist(z2)))
  })
}


getBayesian<-function(w1, w2,wz1,wz2,maxw1,maxw2,ids) {
  corrs<-getZcovariance(wz1[ids],wz2[ids])
  arr=numeric(length(corrs))
  
  for(i in 1:length(corrs)) {
    arr[i]=pnorm(maxw1[i])/pnorm(maxw2[i])
    arr[i]=arr[i]*abs(corrs[i])
  }
  arr
}





plotBarCharts<-function(w1,w2,wz1,wz2){
  
  unique=uniqueBed(w1$peaks,w2$peaks)
  shared=intersectBed(w1$peaks,w2$peaks)
  uniqueid=match(unique$V4,w1$peaks$V4)
  sharedid=match(shared$V4,w1$peaks$V4)  
  maxw1<-w1$getMaxAvgZscore(wz1)
  maxw2<-w2$getMaxAvgZscore(wz2)
  
  
  
  xx=getBayesian(w1,w2,wz1,wz2,maxw1,maxw2,uniqueid)
  yy=getBayesian(w1,w2,wz2,wz1,maxw1,maxw2,sharedid)
  par(mfrow=c(1,2))
  barplot(sort(xx),xlab='conditional probability', ylim=c(0,1))
  title(paste('Unique peaks in',w1$name))
  barplot(sort(yy),xlab='conditional probability', ylim=c(0,1))
  title(paste('Shared peaks in HS959',w1$name))
}


# HSV rainbows
plotMaxAvgZscoreColor<-function(t, w1, w2, wz1,wz2) {
  maxw1<-w1$getMaxAvgZscore(wz1)
  maxw2<-w2$getMaxAvgZscore(wz2)
  id1=w1$peaks$V4
  pvals=getBayesian(w1,w2,wz1,wz2,maxw1,maxw2,id1)
  
  ## Fix data
  pvals[is.nan(pvals)]=0;
  pvals[is.na(pvals)]=0
  pvals[pvals>1]=1
  
  
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  for(i in 1:length(id1))
    points(maxw1[id1][i],maxw2[id1][i],col=hsv(10^(pvals[i]-1)*3/4))
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}


plotMaxAvgZscoreColorUnique<-function(t, w1, w2, wz1,wz2) {
  maxw1<-w1$getMaxAvgZscore(wz1)
  maxw2<-w2$getMaxAvgZscore(wz2)
  unique=uniqueBed(w1$peaks,w2$peaks)
  id2=match(unique$V4,w1$peaks$V4)

  ## Get pvalues and top pvalues
  punique=getBayesian(w1,w2,wz1,wz2,maxw1,maxw2,id2)
  
  ## Fix data
  punique=punique/2.5
  bestp=tail(sort(punique))
  best=match(bestp,punique)
  bestp[bestp>1]=1
  punique[is.nan(punique)]=0;
  punique[is.na(punique)]=0
  punique[punique>1]=1
  
  
  plot(maxw1[id2],maxw2[id2],pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))

  # Color p scores
  for(i in 1:length(id2))
    points(maxw1[id2][i],maxw2[id2][i],col=hsv(punique[i]*3/4))
  
  # Color best p-scores bold
  for(i in 1:length(best)) 
    points(maxw1[id2][best][i],maxw2[id2][best][i],col=hsv(bestp[i]*3/4),cex=2.5,lwd=2)
    
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
  cbind(id2[best],bestp)
}



# HSV rainbows
plotMaxAvgZscoreColorY<-function(t, w1, w2, wz1,wz2) {
  maxw1<-w1$getMaxAvgZscore(wz1)
  maxw2<-w2$getMaxAvgZscore(wz2)
  id1=w1$peaks$V4
  pvals=getBayesian(w1,w2,wz1,wz2,maxw1,maxw2,id1)
  
  ## Fix data
  pvals=pvals/2.5
  pvals[is.nan(pvals)]=0;
  pvals[is.na(pvals)]=0
  pvals[pvals>1]=1
  
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  for(i in 1:length(id1))
    points(maxw1[id1][i],maxw2[id1][i],col=hsv(pvals[i]*3/4))
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}




# HSV rainbows
plotBarChartY<-function(t, wig1, wig2, wz1,wz2) {
  
  unique=uniqueBed(wig1$peaks,wig2$peaks)
  shared=intersectBed(wig1$peaks,wig2$peaks)
  uniqueid=match(unique$V4,wig1$peaks$V4)
  sharedid=match(shared$V4,wig1$peaks$V4)  
  maxw1<-wig1$getMaxAvgZscore(wz1)
  maxw2<-wig2$getMaxAvgZscore(wz2)
  id1=wig1$peaks$V4
  xx=getBayesian(wig1,wig2,wz1,wz2,maxw1,maxw2,uniqueid)
  yy=getBayesian(wig1,wig2,wz1,wz2,maxw1,maxw2,sharedid)
  
  ##pvals=getBayesian(wig1,wig2,wz1,wz2,maxw1,maxw2,id1)
  ##xx=pvals[uniqueid]
  ##yy=pvals[sharedid]
  
  
  ## Fix data

  par(mfrow=c(2,1))
  barplot(sort(xx),ylab='p(x|y)', ylim=c(0,1))
  title('Unique peaks in HS959')
  barplot(sort(yy),ylab='p(x|y)', ylim=c(0,1))
  title('Shared peaks in HS959')
}


plotOverlaps<-function(blist,w1,w2,i,context) {
  bedselect=w2$peaks[blist[i],]
  bedselect[2]=bedselect[2]-context
  bedselect[3]=bedselect[3]+context
  selection1=unlist(w1$Z(bedselect))
  selection2=unlist(w2$Z(bedselect))
  reads1=unlist(w1$getChipReads(bedselect))
  reads2=unlist(w2$getChipReads(bedselect))
  sint=as.numeric(bedselect)
  plotPeakOverlaps(sint,reads1,reads2,selection1,selection2,context)
}

plotPeakOverlaps<-function(bedselect,reads1,reads2,selection1,selection2,context) {
  ## Plot - axis 1
  start=bedselect[2]
  end=bedselect[3]
  l1=length(selection1)
  xx=seq(start-10,end+10,by=10)
  plot(xx,selection1,type='l',xlab='Genome position',ylab='NormDiff score',ylim=c(0,max(selection1)+1),yaxs="i")
  lines(xx,unlist(selection2))
  polygon(c(xx,rev(xx)),c(selection1,numeric(l1)),col='#0000bb77')
  polygon(c(xx,rev(xx)),c(selection2,numeric(l1)),col='#bb000077')
  ## Plot - axis 2
  #par(new=TRUE)
  #l2=length(unlist(reads1))
  #xx=seq(start,end,by=10)
  #plot(xx,reads1,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",lwd=2)
  #lines(xx,reads2,col='red',lwd=2)
  lines(c(start+context,start+context),c(0,1000),col='#00aa0077',lwd=2)
  lines(c(end-context,end-context),c(0,1000),col='#00aa0077',lwd=2)
  title(paste('Peak from', wig1$name, 'vs',wig2$name, 'in chr',bedselect[1]))
  legend('topright',legend=c(wig1$name,wig2$name),fill=c("#0000bb77","#bb000077"))
  #polygon(c(xx,rev(xx)),c(reads1,numeric(l2)),col='#0000bb22')
  #polygon(c(xx,rev(xx)),c(reads2,numeric(l2)),col='#bb000088')
  #axis(4)
  #mtext('ChIP-seq reads',side=4,line=3)
}
#setwd('macs1.4.2')
#e2=loadMacsEnv('S96','HS959')
#setwd('..')
#attach(e2)
#plotTotalReads('Total S96 peak reads vs HS959 synteny MACS 1.4.1 bw=25',wig1,wig2,'pink','green')
#plotTotalReads('Total HS959 peak reads vs S96 synteny MACS 1.4.2 bw=25',wig2,wig1,'lightblue','green')



#####
# e2 is the latest and greatest run....
#plotMaxAvgReads('Max Average reads S96 peaks  w=100', wig1,wig2,'orange','lightgreen')
#plotMaxAvgReads('Max Average reads HS959 peaks  w=100', wig2,wig1,'yellow','lightgreen')

