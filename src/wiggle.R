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
  controlpath=sprintf("%s_MACS_wiggle/control/",name)
  treatpath=sprintf("%s_MACS_wiggle/treat/",name)
  controlname=sprintf("%s_control_afterfiting_",name)
  treatname=sprintf("%s_treat_afterfiting_",name)
  peaks=read.table(sprintf("%s_peaks.bed",name))
  callNextMethod(.Object,...,name=name,
                 controlpath=controlpath,treatpath=treatpath,
                 controlname=controlname,treatname=treatname,
                 peaks=peaks)
})




setGeneric("loadWiggle")

setMethod("loadWiggle","WiggleClass", 
function(this, wigpath, bigwig=TRUE) {
  if(bigwig) {
    files=list.files(path=wigpath,pattern="*.fsa.wig")
    
    #Insertfix or design
    lines<-readLines(files[1])
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
})




setGeneric("loadControlWiggle")
setMethod("loadControlWiggle", "WiggleClass", 
function(this) {
  this@controlwig=this@loadWiggle(this@controlpath,bigwig)
})
setGeneric("loadTreatWiggle")
setMethod("loadTreatWiggle","WiggleClass",
function(this) {
  this@treatwig=this@loadWiggle(this@treatpath,bigwig)
})




##! Delete wig files
#files<-list.files(pattern="*.wig$",recursive=TRUE)
#unlink(files)

##! Unzip all gz files recursively in directory (e.g. for all MACS wiggle files)
#gzip -d -r .






    

