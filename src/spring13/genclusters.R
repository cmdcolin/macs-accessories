## Generate clusters for Cluster and TreeView
# do scaling of data 


out<-apply(normDiffTable[,c(-1,-2)],1,var)
highest<-out>quantile(out,.975)
xout<-cbind(normDiffTable[highest,c(2,1)],normDiffTable[highest,1]+1,var=out[highest])
write.table(xout,row.names=FALSE,col.names=FALSE,sep='\t',file="vartrack.txt",quote=FALSE)



# system.call('bedtools merge -d 20 -i vartrack.txt > varout.txt')




varout<-loadBed('varout.txt',names=c('chr','start','end'))
varout<-varout[(varout[,3]-varout[,2])>50,]
vartab=getPeakScores(varout,normDiffTable,do.rbind=FALSE)

outret<-lapply(vartab,function(x)colSums(x[,c(-1,-2)]))
sumscores<-do.call(rbind,outret)
posret<-lapply(vartab,function(x) c(unique(as.character(x$chr)),min(x[,1]),max(x[,1])))
positions<-do.call(rbind,posret)

strs<-paste0('B',positions[,1],'S',positions[,2],'E',positions[,3])

outret<-lapply(vartab,function(x)colSums(x[,c(-1,-2)]))
sumscores<-do.call(rbind,outret)

sumscoresCol<-as.data.frame(apply(sumscores,2,scale))
write.table(cbind(YORF=strs,sumscoresCol), row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE,file="clusterout-colscale.txt")
sumscoresRow<-t(apply(sumscores,1,scale))
write.table(cbind(YORF=strs,sumscoresRow), row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE,file="clusterout-rowscale.txt")