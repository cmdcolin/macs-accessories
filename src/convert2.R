temp1<-read.table("data/HS959combined_normdiff.txt",skip=1)
temp2<-data.frame(Chr=temp1$V2,temp1$V3,
                         Start=temp1$V3,End=temp1$V3+10,Score=temp1$V5)
write.table(temp2,"out_normdiff.txt",sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)

