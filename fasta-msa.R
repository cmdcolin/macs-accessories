#Fasta binding site sequences
# If this is your first time using the script, install bioconductor using:
#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")


library(ShortRead)


list1=read.table('HS959_only_peaks.bed')
pattern='*.fsa'
fasta=readFasta("S288c-genome/",pattern)


## Make fasta sequences from binding sites
fasta2=c()
idset=c()
for(i in 1:length(list1$V1)) {
  chrnum=as.numeric(substr(list1$V1[i],4,5))
  start = list1$V2[i]
  end=list1$V3[i]
  
  # Get DNA subsequences
  fasta2=c(fasta2,toString(subseq(sread(fasta)[[chrnum]],start,end)));
  # Get ID for binding site
  idset=c(idset,toString(as.character(list1$V4[i])))
}

fastaout=DNAStringSet(fasta2)
ids=BStringSet(idset)
outputfasta=ShortRead(fastaout, ids)
writeFasta(outputfasta, 'fasta-out3.fa')

