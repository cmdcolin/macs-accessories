
#convert all eland_result to bed 
system("for f in *eland_result.txt; do elandresult2bed $f $f.bed; done")

#look for bed files
bedfiles<-list.files(pattern="*.bed")
input.match<-str_match(bedfiles,"(.*)inputu_eland_result.txt.bed")
switchcases<-complete.cases(input.match)
input.match<-input.match[switchcases,]
input.match
chip.match<-bedfiles[!switchcases]
chip.match


#get matches of chip with input files
matches<-lapply(input.match[,2],function(input.sample) {
  match.samples(chip.match,input.sample)
})
names(matches)<-input.match[,2]

match.samples(chip.match,input.match[1,2])
match.samples<-function(chip.samples,input.sample) {
  which(!is.na(str_match(chip.samples,input.sample)))
}




####### load chip files and calculate a top variance without normdiff

bedfiles.reads<-lapply(bedfiles,function(filename) {
  read.BED(filename)
})

chip.list<-lapply(bedfiles.reads,function(elt) bin.chip(elt,1000,zero.filter=F)$chip)

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



