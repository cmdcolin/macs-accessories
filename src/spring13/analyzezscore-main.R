
bed1<-read.table('s96rep1-high_peaks.bed')
bed2<-read.table('s96rep2-high_peaks.bed')
names(bed1)<-c('chr','start','end','name','num')
names(bed2)<-c('chr','start','end','name','num')
bed1$start=as.numeric(bed1$start)
bed1$end=as.numeric(bed1$end)



#############


match=getMatchList(chrnames, macswiggle[[1]]$treat,macswiggle[[2]]$treat)
retable=getScores(match, macswiggle[[1]]$treat,macswiggle[[2]]$treat)
plot(retable$col1,retable$col2,pch='.')
out1<-getPeakScores(bed1,retable)
out2<-getPeakScores(bed2,retable)
points(out1$col1,out1$col2,pch='.',col=2)
points(out2$col1,out2$col2,pch='.',col=3)


###################

cmatch1=getMatchList(chrnames,macswiggle[[1]]$treat,macswiggle[[1]]$control)
cmatch2=getMatchList(chrnames,macswiggle[[2]]$treat,macswiggle[[2]]$control)
retable1=getScores(cmatch1, macswiggle[[1]]$treat,macswiggle[[1]]$control)
retable2=getScores(cmatch2, macswiggle[[2]]$treat,macswiggle[[2]]$control)

#out<-getMatchListMod(chrnames,retable1,retable2)
outg<-getScoresMod(retable1,retable2)

plot(outg$col1,outg$col3,pch='.')




rettab<-getscoresmod(chrnames,macswiggle[[1]]$treat,macswiggle[[1]]$control,macswiggle[[2]]$treat,macswiggle[[2]]$control)

plot(rettab$treat1,rettab$treat2,pch='.')


