# Example preprocessing script.

source('src/wiggle.R')
source('src/zscore.R')

f<-getwd()
setwd('data')

dirs<-list.files(pattern="*_MACS_wiggle")
dirnames<-str_replace_all(dirs,"_MACS_wiggle","")
macswiggle<-lapply(dirnames,function(dirname) {
  printf("Processing %s\n",dirname)
  wig=new("WiggleClass",name=dirname)
  loadWiggles(wig)
  
})

for(i in 1:length(macswiggle)) {
  wig=macswiggle[[i]]
  name=dirnames[[i]]
  scaling<-estimateScalingFactor(wig)
  variance<-estimateVarianceAll(wig,scaling)
  ret<-Zall(wig,scaling,variance)
  
  outtable<-NULL
  for(chr in ret) {
    outtable<-rbind(outtable,chr)
  }
  filename<-sprintf("%s_normdiff.txt",name)
  write.table(outtable,filename,quote=FALSE)
}


#plot(as.vector(ret[[1]][,4]),type="l")
setwd(f)

#--kkkkkkk,,,,kkkj.;/;/;;''''''''''''''''''''



# Convert wiggle file chr names
# 


# ret<-intersectBed(x1,x2)
# ret2<-uniqueBed(x1,x2)
# ret3<-uniqueBed(x2,x1)
# 
# printf("Found %d overlapping peaks\n", nrow(ret))
# printf("Found %d unique peaks in x1\n", nrow(ret2))
# printf("Also found %d unique peaks in x2\n", nrow(ret3))
