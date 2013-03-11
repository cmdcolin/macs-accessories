#accessories

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



