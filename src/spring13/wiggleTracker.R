ntops=as.numeric(rownames(tops))
ntopsTab<-wiggleTable[ntops,c(1,2,2)]
ntopsTab2<-ntopsTab[order(ntopsTab$chr,ntopsTab$pos),]
ntopsTab2[,3]<-ntopsTab2[,3]+20
write.table(ntopsTab2,col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t',file='ntopsTab2.txt')


ret<-read.table('ntopsTab-merge.txt')
b=ntops[1]-50
e=ntops[1]+50
chr=wiggleTable[b,1]
pal=sample(brewer.pal(20,'Dark2'),4)
plot(wiggleTable[b:e,2],wiggleTable[b:e,8],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"))
lines(wiggleTable[b:e,2], wiggleTable[b:e,10],col=pal[2],lwd=2)
lines(wiggleTable[b:e,2],wiggleTable[b:e,4],col=pal[3],lwd=2)
lines(wiggleTable[b:e,2],wiggleTable[b:e,6],col=pal[4],lwd=2)
legend('topright',legend=c('HS959rep2','HS959rep1','S96rep1','S96rep2'),fill=pal)


pal=sample(brewer.pal(20,'Dark2'),4)
plot(normDiffTable[b:e,1],normDiffTable[b:e,8],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"))
lines(normDiffTable[b:e,1], normDiffTable[b:e,10],col=pal[2],lwd=2)
lines(normDiffTable[b:e,1],normDiffTable[b:e,4],col=pal[3],lwd=2)
lines(normDiffTable[b:e,1],normDiffTable[b:e,6],col=pal[4],lwd=2)
legend('topright',legend=c('HS959rep2','HS959rep1','S96rep1','S96rep2'),fill=pal)


r1<-wiggleTable[,4]-wiggleTable[,3]/median(wiggleTable[,3]/wiggleTable[,4])
r2<-wiggleTable[,6]-wiggleTable[,5]/median(wiggleTable[,5]/wiggleTable[,6])
plot((r1+r2)/2,r1-r2,pch='.',xlab='A',ylab='M')
xlist<-(r1+r2)/2
ylist<-(r1-r2)
lm1<-lm(ylist~xlist)
lines(xlist,lm1$fitted)
title('S96 Background subtraction (no normalization)')