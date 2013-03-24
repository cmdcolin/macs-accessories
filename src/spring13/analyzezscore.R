#File created 3/1



# Get chromosomes list
chrnames<-names(macswiggle[[1]]$treat)


## Get matching positions from genome alignments
getmatch<-function(chr,t1,t2) {
    findInterval(t1[[chr]]$V1,t2[[chr]]$V1,all.inside=TRUE)
}

## Apply match score to all chrom
getMatchList<-function(chrlist,t1,t2) {
  sapply(chrlist,getmatch,t1,t2)
}

### Get wig scores for treat only
getscores<-function(matchList,t1,t2) {
  ret<-lapply(names(matchList), function(chr) {
    if(debug)
      printf("processing %s\n", chr);
	  
    currmatch<-matchList[[chr]]
    col1<-t1[[chr]]$V2
    col2<-t2[[chr]]$V2[currmatch]
    data.frame(chr=chr, pos=currmatch, col1=col1,col2=col2)
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}


### Get wig scores for both treat and control
getscoresmod<-function(chrlist,t1,c1,t2,c2) {
  ret<-lapply(chrlist, function(chr) {
    if(debug) {
      printf("processing %s\n", chr);
    }
    temp1=t1[[chr]]
    tempc1=c1[[chr]]
    temp2=t2[[chr]]
    tempc2=c2[[chr]]
    cmatch1=findInterval(temp1$V1,tempc1$V1,all.inside=TRUE)
    cmatch2=findInterval(temp2$V1,tempc2$V1,all.inside=TRUE)
    match=findInterval(temp1$V1,temp2$V1,all.inside=TRUE)
    l1=temp1$V2
    l2=tempc1$V2[cmatch1]
    l3=temp2$V2
    l4=tempc2$V2[cmatch2]
    
    if(debug==TRUE) {
      printf("length %d %d %d %d %d %d %d\n", length(l1),length(l2),length(l3),length(l4),length(match),length(l3[match]),length(l4[match]))
    }
    #currmatch<-matchList[[chr]]
    #col1<-t1[[chr]]$V2
    #col2<-t2[[chr]]$V2[currmatch]
    data.frame(chr=chr, pos=temp1$V1, treat1=l1,control1=l2, treat2=l3[match],control2=l4[match])
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}




getscoresmodappend<-function(wigtab,chrlist,t1,c1,i) {
  ret<-lapply(chrlist, function(chr) {
    if(debug) {
      printf("processing %s\n", chr);
    }
    chrpos=wigtab[wigtab$chr==chr,]$pos
    tempc1=c1[[chr]]
    tempt1=t1[[chr]]
    cmatch1=findInterval(tempt1$V1,tempc1$V1,all.inside=TRUE)
    
    
    match2=findInterval(chrpos,tempt1$V1,all.inside=TRUE)
    l3=tempt1$V2
    l4=tempc1$V2[cmatch1]
    
    if(debug==TRUE) {
      printf("length %d %d %d %d %d\n", length(match2),length(l3),length(l4),length(l3[match2]),length(l4[match2]))
    }
    #currmatch<-matchList[[chr]]
    #col1<-t1[[chr]]$V2
    #col2<-t2[[chr]]$V2[currmatch]
    r1<-data.frame(a=l3[match2],b=l4[match2])
    colnames(r1)<-c(paste0('treat',i),paste0('control',i))
    r1
  })
  
  # from R inferno, Burns (2011)
  x<-do.call('rbind', ret) 
  cbind(wiggles,x)
}





slideMean<-function(x,windowsize=100,slide=1){
  idx1<-seq(1,length(x),by=slide);
  idx1+windowsize->idx2;
  idx2[idx2>(length(x)+1)]<-length(x)+1;
  c(0,cumsum(x))->cx;
  return((cx[idx2]-cx[idx1])/windowsize);
}


##############
# Obsolete below


getPeakScores<-function(bed,scores) {
  chrsplit<-split(scores,factor(scores$chr))
  ret<-apply(bed,1,function(row) {
    s=as.numeric(row[2])
    e=as.numeric(row[3])
    
    chrselect<-strsplit(row[1],'.fsa')[[1]]
    if(debug) {
      printf("Processing %s (%d,%d)\n",chrselect,s,e)
    }
    chrsub<-chrsplit[[chrselect]]
    chrsub[chrsub$pos>s&chrsub$pos<e,]
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}


