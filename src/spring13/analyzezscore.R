



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
# Analyze all chromosomes
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
smoothScatter(scores2$col1,scores2$col2)
plot(scores2$col1,scores2$col2,pch='.')
lm1<-lm(scores2$col2~scores2$col1)
lines(scores2$col1,lm1$fitted)
