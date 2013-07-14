ntops=as.numeric(rownames(x))
ntopsTab<-data.frame(wiggleTable[ntops,c(1,2,2)],order(x$adj.P.Val),order(x$adj.P.Val))
ntopsTab2<-ntopsTab[order(ntopsTab$chr,ntopsTab$pos),]
ntopsTab2[,3]<-ntopsTab2[,3]+20
write.table(ntopsTab2,col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t',file='ntopsTab2.txt')
ntops<-x
ret<-read.table('ntopsTab2.bed')
ret<-ret[(order(ret[,4])),]
ntops<-ret



diff<-wiggleTableScale[x$ID[1],]
b=diff$pos-500
e=diff$pos+60
chr=as.character(diff$chr)

difftable<-difftable[rev(order(difftable$V4)),]
b=difftable[2,2]-300
e=difftable[2,3]+800
chr=as.character(difftable[2,1])
region=ret[ret2[,1]==chr & ret2[,2]>b & ret2[,2]<e,]
region2=ret2[ret2[,1]==chr & ret2[,2]>b & ret2[,2]<e,]
pal=(brewer.pal(6,'BrBG'))
plot(region2[,2],region[,6],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"),ylim=c(0,250))
lines(region2[,2], region[,5],col=pal[2],lwd=2)
lines(region2[,2],region[,4],col=pal[3],lwd=2)
lines(region2[,2],region[,3],col=pal[4],lwd=2)
lines(region2[,2],region[,2],col=pal[5],lwd=2)
lines(region2[,2],region[,1],col=pal[6],lwd=2)
legend('topright',legend=c('HS959rep1','HS959rep2','S96rep1','S96rep2'),fill=pal)
title('Differential peak of S96 vs HS959')

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




wiggleDimeSelection<-wiggleTable[(1:length(dime1$inudge$class))*10,]
out<-wiggleDimeSelection[order(dime1$inudge$fdr),]


ret<-out

b=ret[2,2]-60
e=ret[2,2]+60
chr=ret[2,1]

diff<-wiggleTableScale[x$ID[12],]
myret<-read.table('newreplicatediff-0490.tab.out')
myret<-myret[rev(order(myret$V4)),]

pos<-12
b=myret[pos,2]-2000
e=myret[pos,3]+3000
chr=as.character(myret[pos,1])
myret[pos,]
region=wiggleTableScale[wiggleTableScale[,1]==chr & wiggleTableScale[,2]>b & wiggleTableScale[,2]<e,]
pal=(brewer.pal(4,'BrBG'))
plot(region[,2],region[,3],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"),ylim=c(0,300))
lines(region[,2], region[,4],col=pal[2],lwd=2)
lines(region[,2],region[,6],col=pal[3],lwd=2)
lines(region[,2],region[,7],col=pal[4],lwd=2)
legend('topright',legend=c('S96rep1','S96rep2','HS959rep1','HS959rep2'),fill=pal)
title('Broad differential region S96 vs HS959 (limma)')




#####REVIEW LIMMA Statistics
selection=wiggleTable[x$ID,]
#example region
x[(selection[,1]=='chr15')&(selection[,2]>17400)&(selection[,2]<18000)==TRUE,]





x1<-read.table('../Data/ELAND/S96hs959diff_peaks.bed')
x2<-read.table('../Data/ELAND/Hs959s96diff_peaks.bed')
x1[,1]<-sapply(x1[,1],function(x)strsplit(as.character(x),'.fsa')[[1]])
x2[,1]<-sapply(x2[,1],function(x)strsplit(as.character(x),'.fsa')[[1]])


a1<-1
#ret<-(x1[rev(order(x1[,1])),])
#ret<-x1
b=ret[a1,2]-500

e=ret[a1,3]+100
chr=ret[a1,1]
#b= 43091  
#e=44002
#chr="chr01"
b=ntops[3,2]-100
e=ntops[3,3]+100
chr=ntops[2,1]
chr="chr15"
b=566000
e=568000


myret<-wiggleTable[((1:length(dime1$inudge$fdr))*100)[dime1$inudge$fdr<0.00001],]
b=myret[1,2]-1000
e=myret[1,2]+1000
chr=myret[1,1]
pal=sample(brewer.pal(8,'Dark2'),4)

region=wiggleTable[wiggleTable[,1]==chr & wiggleTable[,2]>b & wiggleTable[,2]<e,]
#region=wiggleTableScale[566000:568000,]
pal=(brewer.pal(4,'RdGy'))
plot(region[,2],region[,8],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"))
lines(region[,2], region[,10],col=pal[2],lwd=2)
lines(region[,2],region[,4],col=pal[3],lwd=2)
lines(region[,2],region[,6],col=pal[4],lwd=2)
legend('topright',legend=c('HS959rep1','HS959rep2','S96rep1','S96rep2'),fill=pal)
title('Significant differential peak S96 vs HS959')
text(b+170,120,paste0("-log10(pvalue)=",ret[a1,5]))



region=voom1$E[wiggleTable[,1]==chr & wiggleTable[,2]>b & wiggleTable[,2]<e,]
region2=wiggleTable[wiggleTable[,1]==chr & wiggleTable[,2]>b & wiggleTable[,2]<e,]
pal=(brewer.pal(4,'RdGy'))
plot(region2[,2],region[,3],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"),ylim=c(0,7))
lines(region2[,2], region[,4],col=pal[2],lwd=2)
lines(region2[,2],region[,1],col=pal[3],lwd=2)
lines(region2[,2],region[,2],col=pal[4],lwd=2)
legend('topright',legend=c('HS959rep1','HS959rep2','S96rep1','S96rep2'),fill=pal)
title('Significant differential peak S96 vs HS959')
text(b+170,120,paste0("-log10(pvalue)=",ret[a1,5]))





b=sortdime[8,2]-2000
e=sortdime[8,3]+500
chr=as.character(sortdime[8,1])
region=wiggleTableScale[wiggleTableScale[,1]==chr & wiggleTableScale[,2]>b & wiggleTableScale[,2]<e,]

pal=(brewer.pal(4,'PiYG'))
plot(region[,2],region[,3],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"),ylim=c(0,180))
lines(region[,2], region[,4],col=pal[2],lwd=2)
lines(region[,2],region[,6],col=pal[3],lwd=2)
lines(region[,2],region[,7],col=pal[4],lwd=2)
legend('topright',legend=c('S96rep1','S96rep2','HS959rep1','HS959rep2'),fill=pal)
title('Significant differential peak S96 vs HS959 (DIME)')
text(b+170,120,paste0("-log10(pvalue)=",ret[a1,5]))





mytable<-read.table('newreplicatediff-0260.tab.out')
sortdime<-mytable[rev(order(mytable$V4)),]
b=sortdime[3,2]-1000
e=sortdime[3,3]+500
chr=as.character(sortdime[3,1])
region=wiggleTableScale[wiggleTableScale[,1]==chr & wiggleTableScale[,2]>b & wiggleTableScale[,2]<e,]

pal=(brewer.pal(4,'RdYlBu'))
plot(region[,2],region[,3],type='l',col=pal[1],lwd=2,ylab='Read score',xlab=paste(chr,"Position"),ylim=c(0,250))
lines(region[,2], region[,4],col=pal[2],lwd=2)
lines(region[,2],region[,6],col=pal[3],lwd=2)
lines(region[,2],region[,7],col=pal[4],lwd=2)
#lines(region[,2],region[,8],col=pal[4],lwd=2)
legend('topright',legend=c('S96rep1','S96rep2','HS959rep1','HS959rep2'),fill=pal)
title('Significant differential peak S96 vs HS959 (limma)')
text(b+170,120,paste0("-log10(pvalue)=",ret[a1,5]))
