# Example preprocessing script.

source('src/wiggle.R')

f<-getwd()
setwd('data')
dataset=list()
for(filename in list.files()) {
  dataset[filename]=new("WiggleClass",name=filename)
}
setwd(f)

#--kkkkkkk,,,,kkkj.;/;/;;''''''''''''''''''''
