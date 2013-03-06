
wiggles<-getscoresmod(chrnames,macswiggle[[3]]$treat,macswiggle[[3]]$control,macswiggle[[1]]$treat,macswiggle[[1]]$control)



bed1<-read.table('s96rep1-high_peaks.bed')
bed2<-read.table('s96rep2-high_peaks.bed')
names(bed1)<-c('chr','start','end','name','num')
names(bed2)<-c('chr','start','end','name','num')
bed1$start=as.numeric(bed1$start)
bed1$end=as.numeric(bed1$end)

plot(wiggles$treat1,wiggles$treat2,pch='.',xlab='S96rep1',ylab='HS959rep1')
title('Raw read counts of S96 vs HS959')