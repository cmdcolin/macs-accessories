

doheatmap<-function(table,granularity=1,Rowv=NA,Colv=NA,scale="none",dist=cor) {
  
  if(granularity!=1) {
    sret<-apply(table[,1:ncol(table)],2,function(x,g){
      slideMean(x,g,g)
    },granularity)
  }
  else {
    sret<-table
  }
  
  heatmap.2(as.matrix(sret),col=redgreen(75), scale=scale,key=TRUE,
            density.info='none',trace='none',Rowv=Rowv,Colv=Colv)
}

doTreeView<-function(table,file=tempfile(),granularity=1) {
  
  if(granularity!=1) {
    lemm<-seq(1,nrow(table),by=granularity)
    coord<-paste0(table$chr[lemm],'P',table$pos[lemm])
    
    temp<-apply(table[,c(-1,-2)],2,function(x,g){
      slideMean(x,g,g)
    },granularity)
    sret<-cbind(coord,temp)
  }
  else {
    sret<-table
  }
  write.table(sret,file=file,row.names=FALSE,quote=FALSE,sep='\t')
}



# Accessory function for plotting large heatmaps
resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
