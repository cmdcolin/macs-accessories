# plot hexbins for NCIS ChIP-seq

chip<-read.BED('data/ELAND-BED/SEG1ChIP_repc.bed')
input<-read.BED('data/ELAND-BED/SEG1input.bed')
ret<-bin.data(chip,input,1000)
ret$input[ret$input>1000]<-1000
ret$chip[ret$chip>8000]<-8000
hexbin1<-hexbin(ret$input,ret$chip,xbins=250)
plot(hexbin1,colramp=rainbow)