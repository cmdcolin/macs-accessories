# Example preprocessing script.
debug=TRUE
source('src/wiggle.R')
source('src/bedfile.R')
source('src/lsos.R')
source('src/spring13/analyzezscore-joinmod.R')


f<-getwd()
setwd('data-bigwigs/')

dirs<-list.files(pattern="*_MACS_wiggle")
dirnames<-str_replace_all(dirs,"_MACS_wiggle","")
macswiggle<-lapply(dirnames,function(dirname) {
  printf("Processing %s\n",dirname)
  wig=new("WiggleClass",name=dirname)
  loadWiggles(wig)
})
names(macswiggle)<-dirnames
setwd(f)

cache('macswiggle')


# Get chromosomes list
chrnames<-names(macswiggle[[1]]$treat)

#Global
nsamples<-length(macswiggle)

# Join wiggle files with matching positions into a table
wiggleTable<-joinWiggleFiles(chrnames, macswiggle)
wiggleTable<-plainNames(wiggleTable)
cache('wiggleTable')





# for(i in 1:length(macswiggle)) {
#   wig=macswiggle[[i]]
#   name=dirnames[i]
#   scaling<-estimateScalingFactor(wig)
#   variance<-estimateVarianceAll(wig,scaling)
#   if(debug) {
#     printf("estimated (global) scaling factor: %f\n", scaling)
#     printf("estimated (global) variance: %f\n", variance)
#   }
#   ret<-Zall(wig,scaling,variance)
#   outtable<-NULL
#   for(chr in ret) {
#     outtable<-rbind(outtable,chr)
#   }
#   filename<-sprintf("%s_normdiff_old.txt",name)
#   write.table(outtable,filename,quote=FALSE)
# }


#plot(as.vector(ret[[1]][,4]),type="l")

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
