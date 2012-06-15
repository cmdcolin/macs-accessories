

#################
#! Constructor
WiggleClass<-function(name) {
  nc=list(
    name=name,
    scaling=1,
    variance=1,
    spacing=10,
    controlpath=paste(name,'/',name,'_MACS_wiggle/control/',sep=''),
    treatpath=paste(name,'/',name,'_MACS_wiggle/treat/',sep=''),
    controlname=paste(name,'_control_afterfiting_',sep=''),
    treatname=paste(name,'_treat_afterfiting_',sep=''),
    wiglist=list(),
    peaks=NULL,
    unique=NULL,
    shared=NULL
    )
  
  ########################
  # Read wiggle files from path into memory and assign filesnames
  nc$loadWiggles=function() {
    loadWiggle<-function(wigpath) {
      files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
      for (filename in files) {
        if(debug)
          print(filename)
        file<-paste(wigpath,filename,sep='')
        wig<-read.table(file, skip=2)
        nc$wiglist[[filename]]=wig
      }
    }
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
      wig=nc$wiglist[[wigfile]]
      corr=findInterval(start:end,wig$V1)
      ret=sum(wig$V2[corr])
      if(debug)
        cat(as.integer(end)-as.integer(start), ' ', ret, '\n')
      ret
    }
    # Apply to chip
    apply(bedfile,1,getTotalReads,nc$treatname)
  } 
  
  # Get avg reads
  nc$getAvgReads = function(bedfile) {
    getAvgReads<-function(x, filepath)
    {  
      chr=x[1]
      start=x[2]
      end=x[3]
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=nc$wiglist[[wigfile]]
      corr=findInterval(start:end,wig$V1)
      b=head(corr,1)
      e=tail(corr,1)
      sum(wig$V2[corr])/(e-b);
    }
    #Apply to treat data
    apply(bedfile,1,getAvgReads,nc$treatname)
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
      wig=nc$wiglist[[wigfile]];
      bstart=start-window
      bend=end
      maxreads=sapply(seq(bstart,end,by=window),function(p){
        b=p
        e=p+window
        corr=findInterval(b:e,wig$V1)
        sum(wig$V2[corr])/window
      })
      ret=max(maxreads)
      if(debug)
        cat('Found max ',ret,' in ',chr,'\n');
      ret
    }
    # Apply to treated data
    apply(bedfile,1,getMaxAvgReads,filepath=nc$treatname,window)
  }

  
  
  
  ####
  # Use all chromosomes for scaling factor
  nc$estimateScalingFactor <- function() {
    sfactor<-function(x) {
      if(debug) {
        cat(x[1], '\t', x[2], '\n')
      }
      treat=nc$wiglist[[x[1]]]
      control=nc$wiglist[[x[2]]]
      corr=findInterval(treat[,1],control[,1])
      control[corr,2]/treat[corr,2]
    }
    files1=list.files(nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(nc$controlpath,pattern="*.fsa.wig.gz")
    ratio_data=apply(cbind(files1,files2),1,sfactor)
    nc$scaling=median(unlist(ratio_data),na.rm=TRUE)
  }
  
  # Sqrt(Aw+Bw/c), w=all
  nc$estimateVarianceAll<-function() {
    files1=list.files(nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(nc$controlpath,pattern="*.fsa.wig.gz")
    getSignal<-function(file){nc$wiglist[[file]]$V2}
    chip_signal=unlist(lapply(files1,getSignal))
    control_signal=unlist(lapply(files2,getSignal))
    #Average signal
    average_chip=mean(chip_signal)
    average_control=mean(control_signal)
    nc$variance=sqrt(average_chip+average_control/nc$scaling^2)
  }
  
  
  ####
  nc$estimateVarianceWindow<-function(xpos, treat,control,ws,corr1,corr2) {
    # find start end positions
    # check for inconsistencies in data
    #
    # checkspacing is correct b window
    #select signals
    #cat(xpos,'\t(', treat$V1[beginning],',',treat$V1[ending],')',nc$matches, '-', nc$unmatches,'\t',
    #    nc$matches1,' ', nc$unmatches1,'(',beginning, ' ', ending, ') (', nc$shifter, ' ', nc$shifter2, '\n')

    b1=findInterval(xpos[1]-ws,corr1)
    e1=findInterval(xpos[1]+ws,corr1)
    b2=findInterval(xpos[2]-ws,corr2)
    e2=findInterval(xpos[2]+ws,corr2)
    chip_signal=treat$V2[b1:e1]
    control_signal=control$V2[b2:e2]
    average_chip=mean(chip_signal)
    average_control=mean(control_signal)
    sqrt(average_chip+average_control/nc$scaling^2)
  }
  
  nc$Zxi<-function(x,treat,control,window,corr1,corr2) {
    pos1=x[1]
    pos2=x[2]
    (treat$V2[pos1]-control$V2[pos2]/nc$scaling)/
      max(nc$estimateVarianceWindow(x,treat,control,window[1],corr1,corr2),
          nc$estimateVarianceWindow(x,treat,control,window[2],corr1,corr2),
          nc$variance)
  }
  
  
  
  # Calculate Z scores over all wiggle files
  nc$Z<-function(bedfile, window=c(10,100)) {
    # Get max average reads over window size
    getZscore<-function(x,f1,f2,window){
      chr=x[1];
      start=as.integer(x[2]);
      end=as.integer(x[3]);
      tf=paste(f1,chr,'.wig.gz',sep='')
      cf=paste(f2,chr,'.wig.gz',sep='')
      treat=nc$wiglist[[tf]]
      control=nc$wiglist[[cf]]
      
      corr1=findInterval(seq(start,end,by=nc$spacing),treat$V1)
      corr2=findInterval(seq(start,end,by=nc$spacing),control$V1)
      if(debug==TRUE)
        cat(chr,'-\t(',start,',',end, ')\n')
      app=cbind(corr1,corr2)
      app=app[-head(app,max(window)),]
      app=app[-tail(app,max(window)),]
      apply(app,1,nc$Zxi,treat,control,window,corr1,corr2)
    }
    
    apply(bedfile,1,getZscore,nc$treatname,nc$controlname, window)
  }
  
  nc$getMaxAvgZscore<-function(wz,ws=10) {
    gmaz<-function(zlist) {
      reads=array()
      for(i in 1:length(zlist)) {
        b=i
        e=i+ws
        reads[i]=mean(zlist[b:e],na.rm=TRUE);
      }
        
      #reads=sapply(1:length(zlist), function(j) {
      #  b=j
      #  e=j+ws
      #  cat(b, ' ', e, '\n')
      #  mean(zlist[b:e],na.rm=TRUE);
      #})
      max(reads)
    }
    sapply(wz,gmaz)
  }
  
  
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}


intersectBed<-function(nc1,nc2) {
  fixer=apply(nc1$peaks,1,function(x){
    chr1=x[1]
    start1=as.integer(x[2])
    end1=as.integer(x[3])
    if(debug)
      cat(chr1,' ',start1, ' ', end1, '\n')
    sub=nc2$peaks[nc2$peaks$V1==chr1,]
    ret=NA
    apply(sub,1,function(y){
      chr2=y[1]
      start2=as.integer(y[2])
      end2=as.integer(y[3])
      
      if(chr1==chr2) {
        #!(AR < BL || BR < AL)
        if(start2 <= end1 && start1 <= end2) ret=x
      }
    })
    ret
  })
  lapply(fixer,na.omit)
}

uniqueBed<-function(nc1,nc2) {
  fixer=apply(nc1$peaks,1,function(x){
    chr1=x[1]
    start1=as.integer(x[2])
    end1=as.integer(x[3])
    if(debug)
      cat(chr1,' ',start1, ' ', end1, '\n')
    sub=nc2$peaks[nc2$peaks$V1==chr1,]
    ret=NA
    apply(sub,1,function(y){
      chr2=y[1]
      start2=as.integer(y[2])
      end2=as.integer(y[3])
      if(chr1==chr2) {
        #AR < BL || BR < AL
        if(end1 <= start2 || end2 <= start1) ret=x
      }
    })
    #fix??
    ret
  })
  fixer
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

