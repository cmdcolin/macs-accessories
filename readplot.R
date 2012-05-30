########################
# Read wiggle files from path into memory and assign filesnames
loadWiggle<-function(wigpath) {
  files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
  for (i in files) {
    file<-paste(wigpath,i,sep='')
    x<-read.table(file, skip=2)
    assign(i,x,inherits=TRUE)
  }
}

#################
#! Constructor
WiggleClass<-function(name) {
  nc=list(
    name=name,
    scaling=1,
    variance=1,
    spacing=10,
    matches=0,
    unmatches=0,
    matches1=0,
    unmatches1=0,
    shifter=0,
    shifter2=0,
    oldname='',
    controlpath=paste(name,'/',name,'_MACS_wiggle/control/',sep=''),
    treatpath=paste(name,'/',name,'_MACS_wiggle/treat/',sep=''),
    controlname=paste(name,'_control_afterfiting_',sep=''),
    treatname=paste(name,'_treat_afterfiting_',sep='')
    )
  
  nc$loadWiggles=function() {
    loadWiggle(nc$treatpath)
    loadWiggle(nc$controlpath)
  }
  #######
  # Get avg reads
  nc$getTotalReads = function(bedfile) {
    getTotalReads<-function(x,filepath){
      chr=x[1]
      start=x[2]
      end=x[3]
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=get(wigfile)
      b=findInterval(start,wig$V1)
      e=findInterval(end,wig$V1)
      sum(wig$V2[b:e])
    }
    # Apply to chip
    apply(bedfile,1,getTotalReads,filepath=nc$treatname)
  } 
  
  # Get avg reads
  nc$getAvgReads = function(bedfile) {
    getAvgReads<-function(x, filepath)
    {  
      chr=x[1]
      start=x[2]
      end=x[3]
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=get(wigfile)
      b=findInterval(start,wig$V1)
      e=findInterval(end,wig$V1)
      sum(wig$V2[b:e])/(e-b);
    }
    #Apply to treat data
    apply(bedfile,1,getAvgReads,filepath=nc$treatname)
  }
  nc$getMaxAvgReads<-function(bedfile, window,inc) {
    # Get max average reads over window size
    getMaxAvgReads<-function(x, filepath, window)
    {
      maxreads=array()
      chr=x[1];
      start=as.integer(x[2]);
      end=as.integer(x[3]);
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=get(wigfile);
      bstart=start-window
      bend=end
      maxreads=sapply(seq(bstart,end,by=window),function(p){
        b=findInterval(p,wig$V1)
        e=findInterval(p+window,wig$V1)
        sum(wig$V2[b:e])/(e-b)
      })
      ret=max(maxreads)
      cat('Found max ',ret,' in ',chr,'\n');
      ret
    }
    # Apply to treated data
    apply(bedfile,1,getMaxAvgReads,filepath=nc$treatname,window)
  }

  
  
  
  ####
  # Use all chromosomes for scaling factor
  nc$estimateScalingFactor <- function() {
    files1=list.files(nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(nc$controlpath,pattern="*.fsa.wig.gz")
    ratio_data=apply(cbind(files1,files2),1,function(x){
      treat=get(x[1])
      control=get(x[2])
      corr=match(treat[,1],control[,1])
      tsig=treat[,2]
      csig=control[,2]
      sapply(corr,function(p){ csig[p]/tsig[p]})
    })
    nc$scaling=median(unlist(ratio_data),na.rm=TRUE)
  }
  
  # Sqrt(Aw+Bw/c), w=all
  nc$estimateVarianceAll<-function() {
    files1=list.files(nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(nc$controlpath,pattern="*.fsa.wig.gz")
    getSignal<-function(file){get(file)$V2}
    chip_signal=unlist(lapply(files1,getSignal))
    control_signal=unlist(lapply(files2,getSignal))
    #Average signal
    average_chip=mean(chip_signal)
    average_control=mean(control_signal)
    nc$variance=sqrt(average_chip+average_control/nc$scaling^2)
  }
  
  
  ####
  nc$estimateVarianceWindow<-function(xpos, treat,control,window) {
    # find start end positions
    # check for inconsistencies in data
    #
    # checkspacing is correct b window
    if(!is.na(treat$V1[xpos-window/nc$spacing]) && 
      treat$V1[xpos]-treat$V1[xpos-window/nc$spacing]==window) {
      beginning=xpos-window/nc$spacing
      nc$matches=nc$matches+1
    } else {
      beginning=findInterval(xpos-window,treat$V1, all.inside=TRUE)
      nc$unmatches=nc$unmatches+1
    }
    if(treat$V1[xpos+window/nc$spacing]-treat$V1[xpos]==window) {
      ending=xpos+window/nc$spacing
      nc$matches=nc$matches+1
    } else {
      ending=findInterval(xpos+window,treat$V1, all.inside=TRUE)
      nc$unmatches=nc$unmatches+1
    }
    if(!is.na(control$V1[beginning-nc$shifter])&&
      treat$V1[beginning]==control$V1[beginning-nc$shifter]) {
      beginning2=beginning
      nc$matches1=nc$matches1+1
    } else {
      beginning2=findInterval(treat$V1[beginning],control$V1,all.inside=TRUE)
      nc$shifter=beginning-beginning2
      nc$unmatches1=nc$unmatches1+1
    } 
    if(!is.na(control$V1[ending+nc$shifter2]) &&
      treat$V1[ending]==control$V1[ending+nc$shifter2]) {
      ending2=ending
      nc$matches1=nc$matches1+1
    } else { 
      ending2=findInterval(treat$V1[ending],control$V1,all.inside=TRUE)
      nc$unmatches1=nc$unmatches1+1
      nc$shifter2=ending2-ending
    }
    #select signals
    #cat(xpos,'\t(', treat$V1[beginning],',',treat$V1[ending],')',nc$matches, '-', nc$unmatches,'\t',
    #    nc$matches1,' ', nc$unmatches1,'(',beginning, ' ', ending, ') (', nc$shifter, ' ', nc$shifter2, '\n')
    chip_signal=treat$V2[beginning:ending]
    control_signal=control$V2[beginning2:ending2]
    average_chip=mean(chip_signal)
    average_control=mean(control_signal)
    sqrt(average_chip+average_control/nc$scaling^2)
  }
  
  nc$Zxi<-function(x,treat,control,window) {
    x2=x
    (treat$V2[x]-control$V2[x]/nc$scaling)/
      max(nc$estimateVarianceWindow(x,treat,control,window[1]),
          #nc$estimateVarianceWindow(x2,treat,control,window[2]),
          nc$variance)
  }
  
  
  
  # Calculate Z scores over all wiggle files
  nc$Z<-function(bedfile, window=c(10)) {
    # Get max average reads over window size
    nc$matches=0
    nc$unmatches=0
    nc$unmatches1=0
    nc$matches1=0
    getZscore<-function(x,f1,f2,window){
      chr=x[1];
      start=as.integer(x[2]);
      end=as.integer(x[3]);
      tf=paste(f1,chr,'.wig.gz',sep='')
      cf=paste(f2,chr,'.wig.gz',sep='')
      treat=get(tf)
      control=get(cf)
      #select interval from wiggle file
      if(as.character(nc$oldname)!=as.character(chr)){
        nc$shifter=0
        nc$shifter2=0
        nc$oldname=chr
        print('here')
      }
      
      b=findInterval(start:end,treat$V1)
      e=findInterval(end,treat$V1)
      cat(chr,'-\t(',start,',',end, ')-\t',b,',',e,'\t',e-b<0,'\n')
      
      V1<-treat$V1[b:e]
      V2<-sapply(b:e,nc$Zxi,treat,control,window)
      cbind(V1,V2)
    }
    
    apply(bedfile,1,getZscore,f1=nc$treatname,f2=nc$controlname,window)
  }
  
  nc$getMaxAvgZscore<-function(Zscore,window=100) {
    acc1<-function(j,window,vpos,vsig) {
      # binary search
      b=findInterval(j,vpos)
      e=findInterval(j+window,vpos)
      mean(vsig[b:e]);
    }
    acc2<-function(x,window,acc1){
      vpos=x[,1]
      vsig=x[,2]
      bstart=head(vpos,1)
      bend=tail(vpos,1)
      cat(bstart,' ',bend, '\n')
      if(bstart==bend)bend=bend+10
      windows=seq(bstart,bend-window/nc$spacing,by=window)
      reads=sapply(windows, acc1,window,vpos,vsig)
      max(reads)
    }
    sapply(Zscore,acc2,window,acc1)
  }
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}





###########


###########
# Get average Z
getMaxAvgZ<-function(bedfile,Zscore,wsize) {
  reads=array();
  chrold=''
  Zchr=NULL
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    # Avoid reloading env variables
    if(as.character(chrold)!=as.character(chr)){
      name=paste(Zscore[['name']],'_treat_afterfiting_',chr,'.wig.gz',sep='')
      name=paste('Z', name)
      
      Zchr=Zscore[[name]]
      print(name)
      chrold=chr
    }    
    
    bstart=start-wsize/10
    bend=end+wsize/10
    maxreads=array()
    # Avoid extract$column in loop
    wigpos=Zchr[,1]
    wigpeak=Zchr[,2]
    
    for(j in seq(bstart,end,by=100)) {
      # binary search
      b=findInterval(j,wigpos)
      e=findInterval(j+wsize,wigpos)
      # Peak window
      peak=wigpeak[b:e];
      maxreads[j]=mean(peak,na.rm=TRUE);
    }
    reads[i]=max(maxreads,na.rm=TRUE);
    print(cat('Found max reads ', reads[i], ' at ', b, ' ', e,'. Used ', (end-bstart)/100, ' windows'))
    
  }
  reads
}



###########
# Get average Z
getAvgZ_WholeGenome<-function(Zscore,wsize) {
  reads=array();
  for(i in Zscore) {
    if(is.null(dim(i))) {next;}
    wigpos=i[,1]
    wigpeak=i[,2]
    bstart=wsize
    bend=length(wigpos)-wsize
    maxreads=array()
    for(j in seq.int(bstart,bend,by=wsize)) {
      # binary search
      b=findInterval(j,wigpos)
      e=findInterval(j+wsize,wigpos)
      
      # Peak window
      peak=wigpeak[b:e];
      
      maxreads[j/wsize]=mean(peak,na.rm=TRUE);
    }
    reads=c(reads,maxreads)
  }
  reads
}




############
# Old get avg Z
getAvgNormDiff<-function(bedfile, path1, path2, scaling_factor, variance) {
  
  reads=array();
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    
    # Reads from S96
    wigfile1=paste(path1,chr,'.wig.gz',sep='');
    wig1=get(wigfile1);
    treat_peak=wig1[wig1$V1>start & wig1$V1<end,];
    wigfile2=paste(path2,chr,'.wig.gz',sep='');
    wig2=get(wigfile2);
    control_peak=wig2[wig2$V1>start & wig2$V1<end,];
  
    treat_peak=as.array(treat_peak$V2)
    control_peak=as.array(control_peak$V2)
    Z=array()
    
    for(j in 1:(end-start)) {
      Z[j]= (treat_peak[j]-control_peak[j]/scaling_factor)/sqrt(variance)
    }
    reads[i]=mean(Z,na.rm=TRUE);
  }
  reads
}






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

