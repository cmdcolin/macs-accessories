# plot hexbins for NCIS ChIP-seq

chip<-read.BED('ELAND-BED/SEG1ChIP_repc.bed')
input<-read.BED('ELAND-BED/SEG1input.bed')
ret<-bin.data(chip,input,1000)
ret$input[ret$input>1000]<-1000
ret$chip[ret$chip>8000]<-8000
hexbin1<-hexbin(ret$input,ret$chip,xbins=250)
plot(hexbin1,colramp=rainbow)


s96rep1<-read.BED('s96rep1.bed')
s96rep2<-read.BED('s96rep2.bed')
s96rep3<-read.BED('s96rep3.bed')
s96input<-read.BED('s96input.bed')
hs959rep1<-read.BED('hs959rep1.bed')
hs959rep2<-read.BED('hs959rep2.bed')
hs959rep3<-read.BED('hs959rep3.bed')
hs959input<-read.BED('hs959input.bed')
s96rep1.bin<-bin.data(s96rep1,s96input,1000,zero.filter=F)$chip
s96rep2.bin<-bin.data(s96rep2,s96input,1000,zero.filter=F)$chip
s96rep3.bin<-bin.data(s96rep3,s96input,1000,zero.filter=F)$chip
s96input.bin<-bin.data(s96rep3,s96input,1000,zero.filter=F)$input
hs959rep1.bin<-bin.data(hs959rep1,hs959input,1000,zero.filter=F)$chip
hs959rep2.bin<-bin.data(hs959rep2,hs959input,1000,zero.filter=F)$chip
hs959rep3.bin<-bin.data(hs959rep3,hs959input,1000,zero.filter=F)$chip
hs959input.bin<-bin.data(hs959rep3,hs959input,1000,zero.filter=F)$input


chip.list<-list(s96rep1.bin,s96rep2.bin,s96rep3.bin,hs959rep1.bin,hs959rep2.bin,hs959rep3.bin)
minlen<-min(sapply(chip.list,length))
lencols<-lapply(chip.list,function(x)x[1:minlen])
chip.frame<-do.call(cbind,lencols)
chip.var<-apply(chip.frame,1,var)