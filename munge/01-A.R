# Example preprocessing script.

source('src/wiggle.R')

f<-getwd()
setwd('data')
dirs<-list.files(pattern="*_MACS_wiggle")
names<-str_replace_all(dirs,"_MACS_wiggle","")
wiggles<-lapply(names,function(name) {
  printf("Processing %s\n",name)
  wig=new("WiggleClass",name=name)
  wig@loadControlWiggle()
  wig@loadTreatWiggle()
  wig
})


setwd(f)

#--kkkkkkk,,,,kkkj.;/;/;;''''''''''''''''''''



# Convert wiggle file chr names
# 
# for(file in c(list.files(recursive=TRUE,pattern="*.wig$"),
#               list.files(recursive=TRUE,pattern="*.bed$"))) 
# {
#   printf("Converting %s\n",file)
#   # Uses sacCer3 chromosome
#   if(saccer==TRUE)
#     convertFileSacCer(file)
#   # Uses S288c chromosome
#   else if(s288c==TRUE)
#     convertFileS288C(file)
# }

# ret<-intersectBed(x1,x2)
# ret2<-uniqueBed(x1,x2)
# ret3<-uniqueBed(x2,x1)
# 
# printf("Found %d overlapping peaks\n", nrow(ret))
# printf("Found %d unique peaks in x1\n", nrow(ret2))
# printf("Also found %d unique peaks in x2\n", nrow(ret3))
