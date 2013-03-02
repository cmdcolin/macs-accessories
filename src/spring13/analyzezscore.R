



r1<-macswiggle[[1]]$treat$chr01
r2<-macswiggle[[2]]$treat$chr01
c1<-macswiggle[[1]]$control$chr01
c2<-macswiggle[[2]]$control$chr01

#NormDiff
match=findInterval(r1$V1,r2$V1)
plot(r1$V2,r2$V2,pch='.')

# Background subtraction
cmatch=findInterval(r1$V1,c1$V1)
cmatch2=findInterval(r2$V1,c2$V1)

p1<-r1$V2-c1[cmatch,2]
p2<-r2$V2-c2[cmatch2,2]
plot(p1,p2[match],pch='.')




# Background subtraction
cmatch=findInterval(r1$V1,c1$V1)
cmatch2=findInterval(r2$V1,c2$V1)

p1<-r1$V2-c1[cmatch,2]
p2<-r2$V2-c2[cmatch2,2]
plot(p1,p2[match],pch='.')




# Background subtraction and scaling
cmatch=findInterval(r1$V1,c1$V1)
cmatch2=findInterval(r2$V1,c2$V1)

m1=median(c1$V2[cmatch]/r1$V2)
m2=median(c2$V2[cmatch2]/r2$V2)
p1<-r1$V2-c1[cmatch,2]/m1
p2<-r2$V2-c2[cmatch2,2]/m2



#Background subtraction and normalization
v1<-mean(r1$V2)+mean(c1$V2[cmatch])/m1^2
v2<-mean(r2$V2)+mean(c2$V2[cmatch2])/m2^2
plot(p1/v1,p2[match]/v2,pch='.',cex=2)



#Minimize local variation
for(i in 1:nrow(r1)) {
  
}






##
# Analyze all chromosomes unmodified only findinterval mapping
chrnames<-names(macswiggle[[1]]$treat)

getmatch<-function(chr,t1,t2) {
  findInterval(t1[[chr]]$V1,t2[[chr]]$V1,all.inside=TRUE)
}

tmatch1<-sapply(chrnames,getmatch,macswiggle[[1]]$treat,macswiggle[[2]]$treat)
tcmatch1<-sapply(chroms, getmatch,macswiggle[[1]]$treat,macswiggle[[1]]$control)
tcmatch2<-sapply(chroms, getmatch,macswiggle[[2]]$treat,macswiggle[[2]]$control)


getscores<-function(matchList,t1,t2) {
  ret<-lapply(names(matchList), function(chr) {
    currmatch<-matchList[[chr]]
    col1<-t1[[chr]]$V2
    col2<-t2[[chr]]$V2[currmatch]
    printf("length %d\t%d\n",length(col1),length(col2))
    data.frame(chr=chr, pos=currmatch, col1=col1,col2=col2)
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}


scores1<-getscores(tcmatch1,macswiggle[[1]]$treat,macswiggle[[1]]$control)
plot(scores1$col1[1:100000],scores1$col2[1:100000],pch='.')
smoothScatter(scores1$col1,scores1$col2)

scores2<-getscores(tmatch1,macswiggle[[1]]$treat,macswiggle[[2]]$treat)
#smoothScatter(scores2$col1,scores2$col2)
plot(scores2$col1,scores2$col2,pch='.')
lm1<-lm(scores2$col2~scores2$col1)
lines(scores2$col1,lm1$fitted)
rm(lm1)




### \
# Perform background subtraction

getmatch<-function(chr,t1,t2) {
  findInterval(t1[[chr]]$V1,t2[[chr]]$V1,all.inside=TRUE)
}

tmatch1<-sapply(chrnames,getmatch,macswiggle[[1]]$treat,macswiggle[[2]]$treat)
tcmatch1<-sapply(chrnames, getmatch,macswiggle[[1]]$treat,macswiggle[[1]]$control)
tcmatch2<-sapply(chrnames, getmatch,macswiggle[[2]]$treat,macswiggle[[2]]$control)

#3333
# analyze medians
getmedian<-function(chrnames,t1,cmatch) {
  sapply(chrnames, function(chr,match) {
    treat1<-t1$treat[[chr]]$V2
    control1<-t1$control[[chr]]$V2[match[[chr]]]
    median(control1/treat1)
  },cmatch)  
}
medl1<-getmedian(chrnames,macswiggle[[1]],tcmatch1)
medl2<-getmedian(chrnames,macswiggle[[2]],tcmatch2)
plot(factor(medl1),factor(medl2))
############



# Do background subtraction
getscores<-function(chrnames,t1,t2) {
  ret<-lapply(chrnames, function(chr) {
    currmatch<-tmatch1[[chr]]
    treat1<-t1$treat[[chr]]$V2
    treat2<-t2$treat[[chr]]$V2
    control1<-t1$control[[chr]]$V2[tcmatch1[[chr]]]
    control2<-t2$control[[chr]]$V2[tcmatch2[[chr]]]
    med1<-median(control1/treat1)
    med2<-median(control2/treat2)
    printf("processing %s median scale %f %f %d\n",chr,med1,med2,i)
    ret1<-treat1-control1/med1
    ret2<-treat2-control2/med2
    data.frame(chr=chr, pos=currmatch, col1=ret1,col2=ret2[currmatch])
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}

scores3<-getscores(chrnames,macswiggle[[1]],macswiggle[[2]])



plot(scores3$col1,scores3$col2,pch='.')
lm1<-lm(scores3$col2~scores3$col1)
lines(scores3$col1,lm1$fitted)
rm(lm1)




###############
#mNormDiff





getnormdiff<-function(chrnames,t1,t2) {
  ret<-lapply(chrnames, function(chr) {
    currmatch<-tmatch1[[chr]]
    treat1<-t1$treat[[chr]]$V2
    treat2<-t2$treat[[chr]]$V2
    control1<-t1$control[[chr]]$V2[tcmatch1[[chr]]]
    control2<-t2$control[[chr]]$V2[tcmatch2[[chr]]]
    med1<-median(control1/treat1)
    med2<-median(control2/treat2)
    printf("processing %s median scale %f %f\n",chr,med1,med2)
    ret1<-treat1-control1/med1
    ret2<-treat2-control2/med2
    nd1<-ret1/sqrt(mean(treat1)+mean(control1)/med1^2)
    nd2<-ret2/sqrt(mean(treat2)+mean(control2)/med2^2)
    data.frame(chr=chr, pos=currmatch, col1=nd1,col2=nd2[currmatch])
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}

scores4<-getnormdiff(chrnames,macswiggle[[1]],macswiggle[[2]])



plot(scores4$col1,scores4$col2,pch='.')
lm1<-lm(scores4$col2~scores4$col1)
lines(scores4$col1,lm1$fitted)
rm(lm1)




bed1<-read.table('s96rep1-new_peaks.bed')
bed2<-read.table('s96rep2-new_peaks.bed')
names(bed1)<-c('chr','start','end','name','num')
names(bed2)<-c('chr','start','end','name','num')
bed1$start=as.numeric(bed1$start)
bed1$end=as.numeric(bed1$end)
apply(bed1,1,function(row) {
  printf("Processing peak %s (%s,%s)\n",row[4],row[2],row[3])
  ret<-strsplit(row[1],'.fsa')[[1]]
  selection1=scores4[scores4$chr==ret&scores4$pos>row[2]&scores4$pos<row[3],]
  points(selection1$col1,selection1$col2,col=2,pch='.')
})
