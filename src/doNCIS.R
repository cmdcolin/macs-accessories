
#convert all eland_result to bed 
system("for f in *eland_result.txt; do elandresult2bed $f $f.bed; done")

#look for bed files
bedfiles<-list.files(pattern="*.bed")
grepinputs<-str_match(bedfiles,"(.*)inputu_eland_result.txt.bed")
switchcases<-complete.cases(grepinputs)
grepinputs<-grepinputs[switchcases,]
input.match<-grepinputs[,1]
chip.prefixes<-grepinputs[,2]
input.match
chip.match<-bedfiles[!switchcases]
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



####### load chip files and calculate a top variance without normdiff


chip.list<-apply(chip.input.matches,1,function(matchrow) {
  cat(paste(input.match[matchrow[1]]," ",chip.match[matchrow[2]],"\n"))
  input.bed<-read.BED(input.match[matchrow[1]])
  chip.bed<-read.BED(chip.match[matchrow[2]])
  bin.data(chip.bed,input.bed,1000,zero.filter=F)$chip
})


####GET NORMDIFF
normdiff.list<-apply(chip.input.matches,1,function(matchrow) {
  cat(paste(input.match[matchrow[1]]," ",chip.match[matchrow[2]],"\n"))
  input.bed<-read.BED(input.match[matchrow[1]])
  chip.bed<-read.BED(chip.match[matchrow[2]])
  nt<-bin.data(chip.bed,input.bed,1000,zero.filter=F)
  getNormDiff(nt$chip,nt$input)
})




#convert to frame
minlen<-min(sapply(chip.list,length))
lencols<-lapply(chip.list,function(x)x[1:minlen])
chip.frame<-do.call(cbind,lencols)
chip.var<-apply(chip.frame,1,var)
threshold<-quantile(chip.var,0.975)
chip.table<-chip.frame[chip.var>threshold,]
write.table(chip.table,file="normdiff-table.txt",sep='\t',quote=F,row.names=F,col.names=F)
system("bedtools merge -d 20 -i normdiff-table.txt > normdiff-merge.txt")
chip.merge.table<-read.table("normdiff-merge.txt",sep="\t")
nrow(chip.merge.table)
nrow(chip.table)

