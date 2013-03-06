


getpeaknormdiff<-function(bed,scores) {
  chrmatch=""
  chrsub=data.frame()
  ret<-apply(bed,1,function(row) {
    s=as.numeric(row[2])
    e=as.numeric(row[3])
    printf("Processing peak %s (%d,%d)\n",row[4],s,e)
    chrselect<-strsplit(row[1],'.fsa')[[1]]
    if(chrmatch!=chrselect) {
      printf("Here %s %s\n",chrmatch,chrselect)
      chrmatch<<-chrselect
      chrsub<<-scores[scores$chr==chrselect,]
    }
    chrsub[chrsub$pos>as.numeric(row[2])&chrsub$pos<as.numeric(row[3]),]
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}



bed1<-read.table('s96rep1-high_peaks.bed')
bed2<-read.table('s96rep2-high_peaks.bed')
names(bed1)<-c('chr','start','end','name','num')
names(bed2)<-c('chr','start','end','name','num')
bed1$start=as.numeric(bed1$start)
bed1$end=as.numeric(bed1$end)

wignormdiff<-data.frame(chr=wiggles$chr,pos=wiggles$pos,nd1=normdiff1,nd2=normdiff2)
retpeaknd<-getpeaknormdiff(bed1,wignormdiff)
retpeakwig<-getpeaknormdiff(bed1,wiggles)
r1<-getpeaknormdiff(bed1,normdiffout)

r2<-apply(r1[,3:6],2,function(list) slideMean(list,100,100))


heatmap(as.matrix(r2),scale='none')


r2<-apply(retpeakwig[,3:10],2,function(list) slideMean(list,100,100))
heatmap(as.matrix(retpeakwig[1:2000,3:10]),scale='none')


testNDs<-log(retpeakwig$treat3/retpeakwig$treat4)

system.time(retlog<-DIME(testNDs))
DIME.plot.fit(retlog)
