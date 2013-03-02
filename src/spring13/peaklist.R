debug=TRUE
source('src/wiggle.R')
source('src/bedfile.R')


par(mfrow=c(1,1))

#load MACS wiggle files
dirs<-list.files(pattern="*_MACS_wiggle")
dirnames<-str_replace_all(dirs,"_MACS_wiggle","")
macswiggle<-lapply(dirnames,function(dirname) {
  printf("Processing %s\n",dirname)
  wig=new("WiggleClass",name=dirname)
  loadWiggles(wig)
  
})



r1=macswiggle[[1]]$treat$chr01
r2=macswiggle[[2]]$treat$chr01
c1=macswiggle[[1]]$control$chr01
c2=macswiggle[[2]]$control$chr01

match=findInterval(r1$V1,r2$V1)
match2=findInterval(r2$V1,r1$V1)

#plot(r1$V2,r2[match,2])
#smoothScatter(r1$V2,r2[match,2])

diffret<-setdiff(r1[match,1],r2$V1)
diffret2<-setdiff(r2[match2,1],r1$V1)

#Plot position of unaligned regions across genomes

old<-par(mfrow=c(1,2))
plot(r1[diffret,1],r2[diffret,1],pch=19)
plot(r2[diffret2,1],r1[diffret2,1],pch=19,col='grey')
par(old)


#Compare peak lists across genomes
rep1<-read.table('S96rep1_peaks.bed')
rep2<-read.table('S96rep2_peaks.bed')
rep2selection<-rep2[rep2$V1=='chr01',]
rep1selection<-rep1[rep1$V1=='chr01',]
ints<-intersectBed(rep1selection,rep2selection)
points(rep2selection$V2-5000,rep2selection$V2+5000,col=2,pch=19)
points(rep1selection$V2+5000,rep1selection$V2-5000,col=4,pch=19)


#Get mappability profiles????
#Compare wiggle scores of unmatched spots in genome
old<-par(mfrow=c(1,2))
y<-r2[diffret,2]
x<-r1[diffret,2]
lm1<-lm(y~x)
plot(x,y,pch=19)
lines(x,lm1$fitted)
y<-r2[diffret2,2]
x<-r1[diffret2,2]
plot(y,x,pch=19,col=2)
lm2<-lm(x~y)
lines(y,lm2$fitted)
par(old)


#wiggle files
rep2selection<-rep2[rep2$V1=='chr01',]
rep1selection<-rep1[rep1$V1=='chr01',]






#SUBTRACT BACKGROUND
plot(r1$V2,r2[match,2],pch='.')
cmatch=findInterval(r1$V1,c1$V1)
cmatch2=findInterval(r2$V1,c2$V1)
plot(r1$V2-c1[cmatch,2],(r2$V2-c2[cmatch2,2])[match],pch='.')



#nORMDIFF
zr1<-r1$V2
zc1<-c1[cmatch,2]
zr2<-r2$V2
zc2<-c2[cmatch2,2]
m1<-median(zc1/zr1,na.rm=TRUE)
m2<-median(zc2/zr2,na.rm=TRUE)
#without normalization
plot(r1$V2-c1[cmatch,2]/m1,(r2$V2-c2[cmatch2,2]/m2)[match],pch='.')
#with normalization

va1<-sqrt(zr1+zc1/(m1^2))
va2<-sqrt(zr2+zc2/(m2^2))

plot((r1$V2-c1[cmatch,2]/m1)/va1,(((r2$V2-c2[cmatch2,2]/m2))/va2)[match],pch='.')





#debug=TRUE
#source('src/wiggle.R')

#r1=macswiggle[[1]]$treat$chr01
#r2=macswiggle[[2]]$treat$chr01
#match=findInterval(r1$V1,r2$V1)
#plot(r1$V2,r2[match,2])
#smoothScatter(r1$V2,r2[match,2])

#diffret<-setdiff(r1[match,1],r2$V1)

#//plot differences across genomes
#plot(r1[diffret,1],r2[diffret,1],pch=19)


#rep1<-read.table('S96rep1_peaks.bed')
#rep1selection<-rep1[rep1$V1=='chr01',]
#rep2<-read.table('S96rep2_peaks.bed')
#rep2selection<-rep2[rep2$V1=='chr01',]
#intersectBed(rep1selection,rep2selection)
#points(rep2selection$V2-5000,rep2selection$V2+5000,col=2,pch=19)
#points(rep1selection$V2+5000,rep1selection$V2-5000,col=4,pch=19)

