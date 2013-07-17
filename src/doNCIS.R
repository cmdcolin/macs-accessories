
#convert all eland_result to bed 
system("for f in *eland_result.txt; do elandresult2bed $f $f.bed; done")
bedfiles<-list.files(pattern="*.bed")
files<-lapply(bedfiles,function(filename) {
  read.BED(filename)
})

input.match<-str_match(names(files),"(.*)input.*")
chip.match<-str_match(names(files),"(.*)ChIP*")

####### load chip files and calculate a top variance without normdiff
chip.list<-lapply(files,function(elt) bin.chip(elt,1000,zero.filter=F)$chip)
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



