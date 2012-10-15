end=as.vector(apply(rep2,1,function(x)
  {
  if(x[9]=='F'){
    as.numeric(as.character(x[8]))+nchar(x[2])
  } else if(x[9]=='R'){
    as.numeric(as.character(x[8]))-nchar(x[2])
  }
}))


out3=data.frame(chr=rep2[,7],start=rep2[,8],end=end,strand=rep2[,9])


write.table(out3,'outfile3.txt')