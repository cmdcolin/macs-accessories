# Example preprocessing script.
debug=TRUE
full_cache=FALSE
source('src/wiggle.R')
source('src/bedfile.R')
source('src/lsos.R')
source('src/spring13/analyzezscore-joinmod.R')



f<-getwd()
setwd('data/')
dirs<-list.files(pattern="*_MACS_wiggle")
dirnames<-str_replace_all(dirs,"_MACS_wiggle","")
macswiggle<-lapply(dirnames,function(dirname) {
  printf("Processing %s\n",dirname)
  wig=new("WiggleClass",name=dirname)
  loadWiggles(wig)
})
setwd(f)

# Get chromosomes list
#Global
nsamples<-length(macswiggle)
chrnames<-names(macswiggle[[1]]$treat)


#names(macswiggle)<-dirnames
if(full_cache) {
  cache('macswiggle')
  cache('nsamples')
  cache('chrnames')
}




# Join wiggle files with matching positions into a table
wiggleTable<-joinWiggleFiles(chrnames, macswiggle)
wiggleTable<-plainNames(wiggleTable)
cache('wiggleTable')
wiggleTable<-wiggleTable[,c(1,2,4,6,8,10,12,14)]
ret<-as.matrix(wiggleTable[,c(3,4,6,7)])
ret2<-wiggleTable[,c(1,2)]
wiggleTableScale<-normalizeMedianValues(ret)
ret<-cbind(ret2,ret)


cmed<-normalizeCustom(ret)
for(i in 1:length(cmed)) {
  ret[,i]=ret[,i]/cmed
}
# 
# normDiffList<-lapply(1:nsamples,function(i) { 
#   pos=i*2
#   str1<-paste0('V',pos-1)
#   str2<-paste0('V',pos)
#   control<-wiggleTable[[str1]]
#   treat<-wiggleTable[[str2]]
#   if(debug) {
#     printf("Processing column%02d (%s, %s)\n",i,str1,str2)
#   }
#   getNormDiff(treat,control)
# })
# 
# 
# 
# # combine rows into table
# normDiffTable<-as.data.frame(do.call(cbind,normDiffList))
# 
# #get pos and chr columns
# normDiffTable<-with(wiggleTable, cbind(pos,chr,normDiffTable))
# 
# 
# #cache('normDiffTable')
# 
# 
# 
# 
# 
# 
# backSubList<-lapply(1:nsamples,function(i) {
#   pos=i*2
#   str1<-paste0('V',pos-1)
#   str2<-paste0('V',pos)
#   #control<-table[[str1]]
#   treat<-wiggleTable[[str2]]
#   
#   m1=median(treat)
#   ret<-treat-m1
#   ret[ret<0]=0
#   ret
# })
# 
# # combine rows into table
# backSubTable<-as.data.frame(do.call(cbind,backSubList))
# 
# #get pos and chr columns
# backSubTable<-with(wiggleTable, cbind(pos,chr,backSubTable))
# rm(backSubList)
# cache('backSubTable')

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
