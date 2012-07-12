# Regex library
library(stringr)

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
    peaks=read.table(paste(name,'/',name,'_peaks.bed',sep='')),
    wiglist=list()
    )
  
  ########################
  # Read wiggle files from path into memory and assign filesnames
  nc$loadWiggles=function(e=environment()) {
    loadWiggle<-function(wigpath,env) {
      files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
      for (filename in files) {
        fn<-paste(wigpath,filename,sep='')
        if(debug)
          cat(filename, '\n')
        if(exists(filename,env))
          wig<-get(filename,env)
        else {
          wig<-read.table(fn, skip=2)
          assign(filename,wig,env)
        }
        nc$wiglist[[filename]]=wig
      }
    }
    
    loadWiggle(nc$treatpath,e)
    loadWiggle(nc$controlpath,e)
  }
  


  nc$getChipReads <- function(bedfile) {
    getChipReadsX<-function(x,filepath){
      chr=x[1]
      start=as.integer(x[2])
      end=as.integer(x[3])
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=nc$wiglist[[wigfile]]
      corr=findInterval(seq(start,end,by=nc$spacing),wig$V1)
      ret=wig$V2[corr]
      if(debug)
        cat(as.integer(end)-as.integer(start), ' ', ret, '\n')
      ret
    }
    # Apply to chip
    apply(bedfile,1,getChipReadsX,nc$treatname)
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
  
  nc$getMaxAvgReads<-function(bedfile, window=100) {
    # Get max average reads over window size
    getMaxAvgReads<-function(x, filepath, window)
    {
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
    sqrt(average_chip+average_control/nc$scaling^2)
  }
  nc$Zxi<-function(x,treat,control,window,corr1,corr2) {
    pos1=x[1]
    pos2=x[2]
    ma=sapply(window,function(w){
      nc$estimateVarianceWindow(x,treat,control,w,corr1,corr2)
    })
    
    ## Normdiff
    A=treat[pos1,2]
    B=control[pos2,2]
    c=nc$scaling
    sigma=max(ma,nc$variance)
    (A-B/c)/sigma
  }
  
  
  
  # Calculate Z scores over all wiggle files
  nc$Z<-function(bedfile, window=c(1,10)) {
    # Get max average reads over window size
    getZscore<-function(x,f1,f2,window){
      chr=x[1];
      start=as.integer(x[2])
      end=as.integer(x[3])
      tf=paste(f1,chr,'.wig.gz',sep='')
      cf=paste(f2,chr,'.wig.gz',sep='')
      treat=nc$wiglist[[tf]]
      control=nc$wiglist[[cf]]
      mw=max(window)
      corr1=findInterval(seq(start-mw,end+mw,by=nc$spacing),treat$V1)
      corr2=findInterval(seq(start-mw,end+mw,by=nc$spacing),control$V1)
      if(debug==TRUE)
        cat(chr,'-\t(',start,',',end, ')\n')
      app=cbind(corr1,corr2)
      apply(app,1,nc$Zxi,treat,control,window,corr1,corr2)
    }
    
    apply(bedfile,1,getZscore,nc$treatname,nc$controlname, window)
  }
  
 
  
  
  
  # Calculate Z scores over all wiggle files
  nc$Zall<-function(window=c(1,10)) {
    files1=list.files(path=nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(path=nc$controlpath,pattern="*.fsa.wig.gz")
    ret=apply(cbind(files1,files2),1,function(f){
      treat<-nc$wiglist[[f[1]]]
      control<-nc$wiglist[[f[2]]]
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
  
  nc$getMaxAvgZscoreAll<-function(wza,ws=100) {
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
  
  
  
  nc$getMaxAvgZscore<-function(wz,ws=10) {
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
  
  
  
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}



intersectBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    chr1=x[1]
    start1=as.integer(x[2])
    end1=as.integer(x[3])
    sublist=nc2[nc2$V1==chr1,]
    ret=apply(sublist,1,function(y){
      chr2=y[1]
      start2=as.integer(y[2])
      end2=as.integer(y[3])
      pn2=y[4]
      (start2<=end1)&&(start1<=end2)
    })
    sum(ret)>0
  })
  
  nc1[selectrows,]
}
uniqueBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    chr1=x[1]
    start1=as.integer(x[2])
    end1=as.integer(x[3])
    pn1=x[4]
    sublist=nc2[nc2$V1==chr1,]
    ret=apply(sublist,1,function(y){
      chr2=y[1]
      start2=as.integer(y[2])
      end2=as.integer(y[3])
      pn2=y[4]
      #AR < BL || BR < AL
      (start2<=end1)&&(start1<=end2)
    })
    sum(ret)==0
  })
  
  nc1[selectrows,]
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

