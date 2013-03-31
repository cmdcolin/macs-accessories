histrv<-hist(wiggleTable[,4],100,plot=FALSE)
histrv$breaks
histrv$density
histrv2<-hist(wiggleTable[,3],100,plot=FALSE)
histrv2$breaks
histrv2$density

plot(histrv$breaks,c(histrv$counts,1),type='l',log="xy")
lines(histrv2$breaks,c(histrv2$counts,1))




out<-apply(normDiffTable[,c(-1,-2)],1,var)
highest<-out>quantile(out,.975)
xout<-cbind(normDiffTable[highest,c(2,1)],normDiffTable[highest,1]+1,var=out[highest])
write.table(xout,row.names=FALSE,col.names=FALSE,sep='\t',file="vartrack.txt",quote=FALSE)
# system.call(') bedtools merge -d 100 -i vartrack.txt > varout.txt')
varout<-loadBed('varout.txt',names=c('chr','start','end'))
varout<-varout[(varout[,3]-varout[,2])>20,]
vartab=getPeakScores(varout,normDiffTable,do.rbind=FALSE)

outret<-lapply(vartab,function(x)colSums(x[,c(-1,-2)]))
sumscores<-do.call(rbind,outret)
posret<-lapply(vartab,function(x) c(unique(as.character(x$chr)),min(x[,1]),max(x[,1])))
positions<-do.call(rbind,posret)

#xout<-apply(vartab[,c(-1,-2)],2,function(x)slideMean(x,5,5))
#iter=seq(1,nrow(xout),by=5)
#strs2<-strs[iter]
strs<-paste0('B',positions[,1],'S',positions[,2],'E',positions[,3])

outret<-lapply(vartab,function(x)colSums(x[,c(-1,-2)]))
sumscores<-do.call(rbind,outret)

#write.table(cbind(strs2,xout), row.names=FALSE,col.names=FALSE,sep='\t',file="clusterout2.txt")
xout2<-xout+abs(min(xout)+1)
write.table(cbind(YORF=strs2,xout2), row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE,file="clusterout2.txt")