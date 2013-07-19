
#convert all eland_result to bed 
system("for f in *eland_result.txt; do elandresult2bed $f $f.bed; done")

#look for bed files
bedfiles<-list.files(path='data',pattern="*.bed")
grep.results<-str_match(bedfiles,"(.*)inputu_eland_result.txt.bed")
switch.cases<-complete.cases(grep.results)
grep.results<-grep.results[switch.cases,]
chip.prefixes<-grep.results[,2] #sample prefixes from grep results
input.match<-bedfiles[switch.cases] #input filenames
chip.match<-bedfiles[!switch.cases] #chip-seq filenames
input.match
chip.match



#use which to find indexof chips that match inputs
match.samples<-function(chip.samples,input.sample) {
  which(!is.na(str_match(chip.samples,input.sample)))
}

#get matches of chip with input files
matches<-lapply(chip.prefixes,function(input.sample) {
  match.samples(chip.match,input.sample)
})




chip.input.matches<-data.frame(a=numeric(0),b=numeric(0))
for(i in 1:length(matches)) {
  for(j in 1:length(matches[[i]])) {
    chip.input.matches<-rbind(chip.input.matches,c(i,matches[[i]][j]))
  }
}
chip.input.matches


####### load chip files and calculate a top variance without normdiff
owd<-setwd('data')
chip.list<-apply(chip.input.matches,1,function(matchrow) {
  cat(paste(input.match[matchrow[1]]," ",chip.match[matchrow[2]],"\n"))
  input.bed<-read.BED(input.match[matchrow[1]])
  chip.bed<-read.BED(chip.match[matchrow[2]])
  bin.data(chip.bed,input.bed,1000,zero.filter=F)$chip
})
setwd(owd)


####GET NORMDIFF
owd<-setwd('data')
normdiff.list<-apply(chip.input.matches,1,function(matchrow) {
  cat(paste(input.match[matchrow[1]]," ",chip.match[matchrow[2]],"\n"))
  input.bed<-read.BED(input.match[matchrow[1]])
  chip.bed<-read.BED(chip.match[matchrow[2]])
  nt<-bin.data(chip.bed,input.bed,1000,zero.filter=F,chr.vec=c('chr01.fsa'))
  getNormDiff(nt$chip,nt$input)
})
setwd(owd)

####


#convert to frame

standardize.length<-function(nlist) {
  minlen<-min(sapply(nlist,length))
  lencols<-lapply(nlist,function(x)x[1:minlen])
  chip.frame<-do.call(cbind,lencols)
}

if(!is.matrix(normdiff.list)) {
  temp<-standardize.length(normdiff.list)
  len<-nrow(temp)
  normdiff.frame<-data.frame(start=seq(0,length.out=len,by=1000),end=seq(1000,length.out=len,by=1000),temp)
}else {
  temp<-normdiff.list
  len<-nrow(temp)
  normdiff.frame<-data.frame(start=seq(0,length.out=len,by=1000),end=seq(1000,length.out=len,by=1000),temp)
}


#get variance
normdiff.var<-apply(normdiff.frame,1,var)
threshold<-quantile(normdiff.var,0.975)
normdiff.select<-normdiff.frame[normdiff.var>threshold,]

options(scipen=999)
write.table(normdiff.select,file="normdiff-table.txt",sep='\t',quote=F,row.names=F,col.names=F)
options(scipen=0)
system("bedtools merge -d 20 -i normdiff-table.txt > normdiff-merge.txt")
chip.merge.table<-read.table("normdiff-merge.txt",sep="\t")
nrow(chip.merge.table)
nrow(chip.table)

