# WiggleClass is a S4 object that loads Wiggle tracks

setClass("WiggleClass",
         representation(
           name="character", 
           wigpath="character",
           peaks="data.frame",
           wiggledata="list")
)


setMethod("initialize","WiggleClass",function(.Object,...,name) {
  printf("Loading %s\n", name)
  wigpath=sprintf("%s_MACS_wiggle/",name)
  peaks=read.table(sprintf("%s_peaks.bed",name))
  callNextMethod(.Object,...,
                 name=name,
                 wigpath=wigpath, 
                 peaks=peaks)
})



setGeneric("loadBigWig", function(this){
  standardGeneric("loadBigWig")})

setMethod("loadBigWig","WiggleClass", 
function(this) 
{
  files=list.files(path=this@wigpath,pattern="*.wig",recursive=TRUE)
  ret<-lapply(files,function(file) {
    #Insertfix or design
    print(file)
    path<-sprintf("%s/%s",this@wigpath,file)
    lines<-readLines(path)
    ntable<-grep("variableStep",lines)
    ltable<-length(ntable)-1
    
    # Extract each table in bigwig file
    rep<-lapply(1:ltable, function(i) {
      varline<-lines[ntable[i]]
      chrom<-unlist(str_extract_all(varline,"[Cc]hr[0-9a-z]*"))[2]
      if(debug) {
        printf("Processing %s\n", chrom)
      }
      
      begin=ntable[i]+1
      end=ntable[i+1]-1
      data<-paste(lines[begin:end],collapse='\n')
      con<-textConnection(data)
      tab5rows <- read.table(con, nrows = 5)
      classes <- sapply(tab5rows, class)
      chr <- read.table(con, colClasses = classes,nrows=length(lines))
      close(con)
      attr(chr, 'name') <- chrom  # save the nfame
      chr
    })
    names(rep)<-sapply(rep,attr,'name')
    rep
  })
  names(ret)<-c("control","treat")
  ret
})


setGeneric("loadWiggles", function(this, bigwig=TRUE){
  standardGeneric("loadWiggles")})

# Load Wiggles
setMethod("loadWiggles","WiggleClass", 
function(this, bigwig=TRUE) 
{
  if(bigwig) {
    this@wiggledata=loadBigWig(this)
  }
  else {
    files=list.files(pattern="*.wig",recursive=TRUE)
    lapply(files,function(filename) {
      if(debug)
        cat(filename, '\n')
      wig<-read.table(filename, skip=2, header=TRUE)
    })
  }
  this@wiggledata
})




# Set wiggle table names as abbreviated sample names
prettyNames<-function(wiggleTable) {
  
  fixNames<-sapply(names(macswiggle),function(x) strsplit(x,'-new')[[1]])
  prettyNames<-paste0(sort(rep(fixNames,2)),'-',rep(c('c','t'),nsamples))
  names(wiggleTable)<-c('chr','pos',prettyNames)
}

## Set wiggle table names as chr, pos, V1-Vn
plainNames<-function(wiggleTable) {
  
  fixNames<-c('chr','pos',paste0(rep('V',nsamples),1:(nsamples*2)))
  names(wiggleTable)<-fixNames
  wiggleTable
}
#setGeneric("loadWiggles", function(this, bigwig=TRUE){
#  standardGeneric("loadWiggles")})
#setMethod("loadWiggles", "WiggleClass", 
#function(this, bigwig=TRUE) {
#  this@controlwig=loadWiggle(this, this@controlpath, bigwig)
#  this@treatwig=loadWiggle(this, this@treatpath, bigwig)
#})




##! Delete wig files
#files<-list.files(pattern="*.wig$",recursive=TRUE)
#unlink(files)

##! Unzip all gz files recursively in directory (e.g. for all MACS wiggle files)
#gzip -d -r .






    

