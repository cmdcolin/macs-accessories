for(i in 1:500) {
  ret<-topTags(et,n=100+10*i)
  sel=wiggleTable[rownames(ret),]
  sel2<-data.frame(sel[,1],sel[,2],sel[,2]+10,sel[,c(-1,-2)])
  sel3<-sel2[order(sel2[,1],sel[,2]),]
  
  write.table(sel3,quote=F,row.names=F,col.names=F,file=paste0("../wiggleTabSel",i,".bed"),sep='\t')
}
