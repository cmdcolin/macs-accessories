# Example preprocessing script.

source('src/wiggle.R')

f<-getwd()
setwd('data')
dirs<-list.files(pattern="*_MACS_wiggle")
store<-list()
x<-list()
macswiggle<-lapply(dirs,function(dir) {
  name<-str_replace(dir,"_MACS_wiggle","")
  printf("Processing %s\n",name)
  wig=new("WiggleClass",name=name)
  loadWiggles(wig)
  x<-wig
  scaling<-estimateScalingFactor(wig)
  variance<-estimateVarianceAll(wig,scaling)
  ret<-Zall(wig,scaling,variance)
  write.table(ret,sprintf("%s_normdiff.txt",name))
})

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
