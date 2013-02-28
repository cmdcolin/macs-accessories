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

