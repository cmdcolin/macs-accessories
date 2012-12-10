nfiles=list.files(pattern="*_normdiff.txt")
bfiles=list.files(pattern="*_peaks.bed")
ret=lapply(nfiles,function(x) read.table(x))
ret2=lapply(bfiles,function(x) read.table(x))
hrep1=ret[[1]]
hrep1_peaks=ret2[[1]]
chr1=hrep1[hrep1$V1=='chr01',]
chr1_peaks=hrep1_peaks[hrep1_peaks$V1=='chr01',]
plot(chr1$V2,chr1$V4,type='l')
ret=apply(chr1_peaks,1,function(x){
  rect(x[2],-100,x[3],100,col=rgb(1,0,0,0.5))
})



require(R.utils)
qqnd<-function(normdiff) {
  clt=sapply(1:5000,function(x){mean(sample(normdiff,50))})
  nd.mu=mean(clt)
  nd.sd=sd(clt)
  qqnorm(clt)
  qqline(rnorm(1000,mean=nd.mu,sd=nd.sd))
}
qqnd(normdiff=hrep1$V4)

densitynd<-function(normdiff) {
  clt=sapply(1:10000,function(x){mean(sample(normdiff,50))})
  nd.mu=mean(clt)
  nd.sd=sd(clt)
  plot(function(x){dnorm(x,mean=nd.mu,sd=nd.sd)},-3,3)
  lines(density(clt))
  printf("%f\t%f",nd.mu,nd.sd)
  #lines(density(rnorm(1000,nd.mu,nd.sd)))
}
densitynd(normdiff=hrep1$V4)



#################
#get all chr names
#loop from seq 0 to maxchr length by window size
#match intervals in samples for samples in location
#average normdiff scores
#plot pair on plot
##
nd1=ret[[1]]
nd2=ret[[2]]
str=c('chr','start','end','score')
names(nd1)<-str
names(nd2)<-str
plot(0,0,xlim=c(-3,15),ylim=c(-3,15))
chrs=levels(nd1$chr)
results=lapply(chrs,function(name) {
  exchr1=nd1[nd1$chr==name,]
  exchr2=nd2[nd2$chr==name,]
  
  intervals=findInterval(as.integer(exchr1$start),as.integer(exchr2$start),all.inside=TRUE)
  len=length(intervals)
  points=seq(1,len-10,by=10)
  x=sapply(intervals[points],function(i) {
    pos=findInterval(i,intervals)
    sel=intervals[pos:(pos+10)]
    
    m1=mean(exchr1[sel,]$score)
    m2=mean(exchr2[sel,]$score)
    points(m1,m2)
    c(m1,m2,pos)
  }) 
})

# gen all points
chr1=results[[1]]
plot(chr1[1,],chr1[2,],pch='.')
for(i in 1:16) {
  chr=results[[i]]
  points(chr[1,],chr[2,],pch='.')
}



# find intersections
bed1=ret2[[1]]
bed2=ret2[[2]]
str=c('chr','start','end','name','score')
names(bed1)<-str
names(bed2)<-str

for(i in 1:16) {
  name=chrs[i]
  exchr1=nd1[nd1$chr==name,]
  exchr2=nd2[nd2$chr==name,]
  exbed1=bed1[bed1$chr==name,]
  exbed2=bed2[bed2$chr==name,] 
  chr=results[[i]]
  list1=unlist(apply(exbed1,1,function(x){seq(as.integer(x[2]),as.integer(x[3]),10)}))
  list2=unlist(apply(exbed2,1,function(x){seq(as.integer(x[2]),as.integer(x[3]),10)}))
  intervals1=findInterval(list1,chr[3,],all.inside=TRUE)
  intervals2=findInterval(list2,chr[3,],all.inside=TRUE)
  print(str(exbed1$start))
  print(str(exbed2$start))
  print(str(chr[3,intervals1]))
  print(str(chr[3,intervals2]))
  print(str(intervals1))
  print(str(intervals2))
  printf("\n\n\n")
  
  
  points(chr[1,intervals1],chr[2,intervals1],pch=20,col=rgb(1,0,0,0.5))
  points(chr[1,intervals2],chr[2,intervals2],pch=20,col=rgb(0,0,1,0.5))
}

