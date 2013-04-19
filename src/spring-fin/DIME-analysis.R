

temp1<-(1:length(dime1$best$class))*10
temp2<-wiggleTable[temp1[dime1$best$class==1],]
write.table(data.frame(temp2[,1],temp2[,2],temp2[,2]+10),row.names=FALSE,col.names=FALSE,quote=FALSE,file="ndimeTab.bed",sep='\t')
ndimeTab<-read.table('ndimeTab-merge.bed')
# topTable limma merged
ntopsTab<-read.table('ntopsTab2.bed')