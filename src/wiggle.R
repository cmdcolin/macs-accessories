# Use regex library
library(stringr)
library(R.utils)


setClass("WiggleClass",
         representation(
           name="character", 
           controlpath="character",
           treatpath="character",
           controlname="character",
           treatname="character",
           peaks="data.frame",
           controlwig="list",
           treatwig="list")
)


setMethod("initialize","WiggleClass",function(.Object,...,name) {
  printf("Loading %s\n", name)
  controlpath=sprintf("%s_MACS_wiggle/control/",name,name)
  treatpath=sprintf("%s_MACS_wiggle/treat/",name,name)
  controlname=sprintf("%s_control_afterfiting_",name,name)
  treatname=sprintf("%s_treat_afterfiting_",name,name)
  peaks=read.table(sprintf("%s_peaks.bed",name,name))
  callNextMethod(.Object,...,name=name,
                 controlpath=controlpath,treatpath=treatpath,
                 controlname=controlname,treatname=treatname,
                 peaks=peaks)
})

WiggleClass.loadWiggle<-function(.Object, wigpath, bigwig=TRUE) {
  if(bigwig) {
    lines<-readLines(file)
    ntable<-grep("variableStep",lines)
    
    lapply(1:length(n), function(i) {
      
      if(debug)
        printf("%s\n", ntable[i])
      
      begin=ntable[i]+1
      end=ntable[i+1]-1
      con<-textConnection(lines[begin:end])
      chr<-read.table(con, skip=1, header=TRUE)
      close(con)
      
      chr
    })
  }
  else {
    files=list.files(path=wigpath,pattern="*.fsa.wig")
    ret<-lapply(files,function(x){
      filename<-paste(wigpath,filename,sep='')
      if(debug)
        cat(filename, '\n')
      wig<-read.table(filename, skip=2, header=TRUE)
    })
    ret
  }
}

setMethod("loadWiggle","WiggleClass", WiggleClass.loadWiggle)



# Where are my ROMAN NUMERAL FUNCTIONS??
# Can't upload anything to RESEARCH COMPUTING??
# Evaluate MANORM pairwise on replicates
# Concatenate replicate data into one file
# Compare results
# Fix many roman numeral sequences
#Figure out correct regex

romanNum<-function(str) {
  if(str=="01") "I"
  else if(str=="02") "II"
  else if(str=="03") "III"
  else if(str=="04") "IV"
  else if(str=="05") "V"
  else if(str=="06") "VI"
  else if(str=="07") "VII"
  else if(str=="08") "VIII"
  else if(str=="09") "IX"
  else if(str=="10") "X"
  else if(str=="11") "XI"
  else if(str=="12") "XII"
  else if(str=="13") "XIII"
  else if(str=="14") "XIV"
  else if(str=="15") "XV"
  else if(str=="16") "XVI"
  else if(str=="mt") "M"
  else "Error"
}



# Uses roman numerals
convertFileSacCer3<-function(filename) {
  
  filetext<-readLines(filename)
  filetext <- paste(filetext,collapse="\n")
  
  for(i in 1:16) {
    str<-sprintf("%02d",i)
    filetext<-str_replace_all(filetext, 
      sprintf("chr%s.fsa",str), sprintf("chr%s",romanNum(str)))
    printf("Finished chr%02d\n", i)
  }
  filetext<-str_replace_all(filetext,  "chrmt.fsa", "chrM")
  printf("Finished chrmt\n", i)
  writeLines(filetext,filename)
  
}

convertFileS288C<-function(file) {
  
  filetext<-readLines(file)
  filetext <- paste(filetext,collapse="\n")
  
  for(i in 1:16) {
    str<-sprintf("%02d",i)
    filetext<-str_replace_all(filetext, 
      sprintf("chr%s.fsa",str), sprintf("Chr%s",str))
    if(debug==TRUE)
      printf("Finished chr%02d\n", i)
  }
  filetext<-str_replace_all(filetext,  "chrmt.fsa", "Chrmt")
  if(debug==TRUE)
    printf("Finished chrmt\n", i)
  
  writeLines(filetext,file)
  
}



##! Delete wig files
#files<-list.files(pattern="*.wig$",recursive=TRUE)
#unlink(files)

##! Unzip all gz files recursively in directory (e.g. for all MACS wiggle files)
#gzip -d -r .



loadBed<-function(path) {
  bedfile<-read.table(path)
  
  # rename columns
  names(bedfile)<-c("chromosome","start","end","name","value")
  
  #make numeric columns
  transform(bedfile, {
    start = as.numeric(as.character(start))
    end = as.numeric(as.character(end))
  })
}



intersectBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    
    sublist=nc2[nc2$chromosome==x['chromosome'],]
    ret=apply(sublist,1,function(y){
      #  overlap AR < BL || BR < AL
      (y['start']<=x['end'])&&(x['start']<=y['end'])
    })
    
    ##! Get overlap peaks where intersect>0
    sum(ret)>0
  })
  nc1[selectrows,]
}



intersectBedLimma<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    
    sublist=nc2[nc2$chromosome==x['chromosome'],]
    ret=apply(sublist,1,function(y){
      #  overlap AR < BL || BR < AL
      (y['start']<=x['end'])&&(x['start']<=y['end'])
    })
    
    ##! Get overlap peaks where intersect>0
    sum(ret)>0
  })
  selectrows
}


uniqueBed<-function(nc1,nc2) {
  selectrows=apply(nc1,1,function(x){
    
    sublist=nc2[nc2$chromosome==x['chromosome'],]
    ret=apply(sublist,1,function(y){
      # Linear overlap AR < BL || BR < AL
      (y['start']<=x['end'])&&(x['start']<=y['end'])
    })
    
    ##! Get unique peaks with no overlap
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

    
    
###! Unused

    #extended regex
    #str_replace_all(text, "(chr)([0-9a-z]{2})(.fsa)", sprintf("\\1%s","\\1",getRoman("\\2")))
convertFileOld<-function(filename) {
  con<-file(filename)
  open(con)
  lines<-readLines(con,-1)
  outlines<-sapply(lines, function(cline) {
    x<-str_match(cline,"(.*chr)([0-9a-z]{2})(.fsa)(.*)")
    if(sum(is.na(x))==0) {
      sprintf("%s%s%s",x[2],romanNum(x[3]),x[5])
    }
    else
      cline
  })
  writeLines(outlines,filename)
  close(con)
}
