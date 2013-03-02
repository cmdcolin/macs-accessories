chrnames<-names(macswiggle[[1]]$treat
                

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