bin.size<-20

#for f in normdiff*; do bedtools merge -d 200 -i $f -scores sum > $window.size<-1000
#convert all eland_result to bed 
#system("for f in *eland_result.txt; do elandresult2bed $f $f.bed; done")

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
# owd<-setwd('data')
# chip.list<-apply(chip.input.matches,1,function(matchrow) {
#   cat(paste(input.match[matchrow[1]]," ",chip.match[matchrow[2]],"\n"))
#   input.bed<-read.BED(input.match[matchrow[1]])
#   chip.bed<-read.BED(chip.match[matchrow[2]])
#   bin.data(chip.bed,input.bed,1000,zero.filter=F)$chip
# })
# setwd(owd)


####GET NORMDIFF

#load bed files
owd<-setwd('data')
bedfiles<-apply(chip.input.matches,1,function(matchrow) {
  cat(paste(input.match[matchrow[1]]," ",chip.match[matchrow[2]],"\n"))
  input.bed<-read.BED(input.match[matchrow[1]])
  chip.bed<-read.BED(chip.match[matchrow[2]])
  list(chip.bed,input.bed)
})
cache('bedfiles')
#create binned data by chromosome
chippies<-lapply(bedfiles,function(row) {
  nt<-bin.data(row[[1]],row[[2]],bin.size,zero.filter=F,by.chr=TRUE)
})
cache('chippies')


#calculate normdiff
normdiff.list<-lapply(chippies,function(nt){
  chrnames<-names(nt$chip)#assume that noth input,chip have same names
  ret<-lapply(chrnames,function(chr) {
    getNormDiff(nt$chip[[chr]],nt$input[[chr]])
  })
  names(ret)<-names(nt$chip)
  ret
})
names(normdiff.list)<-chip.prefixes[chip.input.matches[,1]]
setwd(owd)
#use filenames
#names(normdiff.list)<-chip.match[chip.input.matches[,2]]
#use prefixes

####


#convert to frame
standardize.length<-function(nlist,chrnames) {
  print(names(nlist))
  ret<-lapply(chrnames,function(chrname) {
    
    min.length<-min(sapply(nlist,function(chipseqbins) {
      length(chipseqbins[[chrname]])
    }))
    columns<-lapply(nlist,function(chrdata) {
      chrdata[[chrname]][1:min.length]
    })
    
    data.frame(chr=chrname,
               start=seq(0,length.out=min.length,by=bin.size),
               end=seq(bin.size,length.out=min.length,by=bin.size),
               do.call(cbind,columns))
  })
  
  do.call(rbind,ret)
}

normdiff.frame<-standardize.length(normdiff.list,names(normdiff.list[[1]]))
print(nrow(normdiff.frame))




#get variance

normdiff.var<-apply(as.matrix(normdiff.frame[,c(-1,-2,-3)]),1,var)
#normdiff.var=numeric(nrow(normdiff.frame))
#for(i in 1:nrow(normdiff.frame)){normdiff.var[i]=var(as.numeric(normdiff.frame[i,c(-1,-2,-3)]))}

#threshold variances
threshold<-quantile(normdiff.var,0.975)
normdiff.select<-normdiff.frame[normdiff.var>threshold,]
print(nrow(normdiff.select))



#output as bed file
options(scipen=999) #remove scientific notation


for(i in 4:ncol(normdiff.select)) {
  filename<-sprintf("normdiff-table%02d.txt",i)
  cat(paste("writing",filename,"\n"))
  write.table(normdiff.select[,c(1,2,3,i,i)],file=filename,sep='\t',quote=F,row.names=F,col.names=F)
}
options(scipen=0)

#merge files:
#for f in normdiff*.txt; do bedtools merge -d 200 -i $f -scores sum > $f.merge; done
#

bedfiles<-list.files(pattern="*.merge")
normdiff.merge<-read.table(bedfiles[1])
for(i in 2:length(bedfiles)) {
  normdiff.merge<-cbind(normdiff.merge,read.table(bedfiles[i])[,4])
}
names(normdiff.merge)<-c('chr','start','end',chip.prefixes[chip.input.matches[,1]])


#merge
heatmap.2(as.matrix(normdiff.merge[,c(-1,-2,-3)]),Rowv=NA,trace="none",scale="none")
boxplot(as.matrix(normdiff.merge[,c(-1,-2,-3)]),boxwex=0.6, notch=T, outline=FALSE,las=2)

read.depth.table<-read.table('readdepth.txt')
ret<-str_match(read.depth.table[,4],"(.*)input")
cbind(chip.prefixes[chip.input.matches[,1]],read.depth.table[!complete.cases(ret),1][1:17])
cbind(chip.match[chip.input.matches[,2]],read.depth.table[!complete.cases(ret),1][1:17])

     
system("bedtools merge -d 20 -i normdiff-table.txt > normdiff-merge.txt")
chip.merge.table<-read.table("normdiff-merge.txt",sep="\t")
nrow(chip.merge.table)
nrow(chip.table)

