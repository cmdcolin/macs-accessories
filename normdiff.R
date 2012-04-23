# April 15th 2012
# 
#
# Change log 4/21/12 use median in estimate_scaling_factor
#
#
#

estimate_scaling_factor <- function(treat, control) {
  ratio_data=as.array(t(control$V2/treat$V2))
  ratio_data_median=apply(ratio_data, 1, median)
  ratio_data_median
}

# Sqrt(Aw+Bw/c), w=all
estimate_variance_all<-function(treat,control,scaling_factor) {
  chip_signal=as.array(t(treat$V2))
  control_signal=as.array(t(control$V2))
  average_chip=apply(chip_signal,1,mean)
  average_control=apply(control_signal,1,mean)
  sqrt(average_chip+average_control/scaling_factor^2)
}

# Sqrt(Aw+Bw/c), w=1
estimate_variance_one<-function(treat,control,scaling_factor,pos) {
  chip_signal=as.array(t(treat$V2[(pos-1):(pos+1)]))
  control_signal=as.array(t(control$V2[(pos-1):(pos+1)]))
  average_chip=apply(chip_signal,1,mean)
  average_control=apply(control_signal,1,mean)
  sqrt(average_chip+average_control/scaling_factor^2)
}

# Sqrt(Aw+Bw/c), w=10
estimate_variance_ten<-function(treat,control,scaling_factor,pos) {
  chip_signal=as.array(t(treat$V2[(pos-10):(pos+10)]))
  control_signal=as.array(t(control$V2[(pos-10):(pos+10)]))
  average_chip=apply(chip_signal,1,mean)
  average_control=apply(control_signal,1,mean)
  sqrt(average_chip+average_control/scaling_factor^2)
}



max_estimate_variance<-function(treat,control,scaling_factor,pos,varianceall) {
  max(estimate_variance_one(treat,control,scaling_factor,pos),
      estimate_variance_ten(treat,control,scaling_factor,pos),
      varianceall)
}


Zxi<-function(treat,control,scaling_factor,pos,varianceall) {
  (treat$V2[pos]-control$V2[pos]/scaling_factor)/
    max_estimate_variance(treat,control,scaling_factor,pos,varianceall)
}


# Read data files
setwd('C:/Documents and Settings/Colin Diesh/Desktop/Workspace/')
treat=read.table('S96/S96_MACS_wiggle/treat/S96_treat_afterfiting_chr01.fsa.wig.gz', skip=2)
control=read.table('S96/S96_MACS_wiggle/control/S96_control_afterfiting_chr01.fsa.wig.gz', skip=2)
treat=treat[-length(treat$V2),]

# Prepare variances
scaling_factor=estimate_scaling_factor(treat,control)
varianceall=estimate_variance_all(treat,control,scaling_factor)


##  Get standard deviations
start=10
end=length(treat$V2)-10
Z=lapply(start:end, function(x){Zxi(treat,control,scaling_factor,x,varianceall)})




plot(1:1000,Z[1:1000],type='l', ylab='Z-score', xlab='Genome position')
plot(1000:2000,Z[1000:2000],type='l', ylab='Z-score', xlab='Genome position')
plot(2000:3000,Z[2000:3000],type='l', ylab='Z-score', xlab='Genome position')
plot(3000:4000,Z[3000:4000],type='l', ylab='Z-score', xlab='Genome position')
plot((4000:5000)*10,Z[4000:5000],type='l', ylab='Z-score', xlab='Genome position')






######################

#NOTE length of files different. Handle somehow...
# MACS probably doesn't report locations where there are 0 reads
# Read data files
# Okay, use commas for each dimension of table
treat=read.table('HS959/HS959_MACS_wiggle/treat/HS959_treat_afterfiting_chr01.fsa.wig.gz', skip=2)
control=read.table('HS959/HS959_MACS_wiggle/control/HS959_control_afterfiting_chr01.fsa.wig.gz', skip=2)
treat=treat[-(length(treat$V2):(length(control$V2)+1)),]
length(treat$V2)
length(control$V2)
# Prepare variances
scaling_factor=estimate_scaling_factor(treat,control)
varianceall=estimate_variance_all(treat,control,scaling_factor)


##  Get standard deviations
start=10
end=length(treat$V2)-10
Z=lapply(start:end, function(x){Zxi(treat,control,scaling_factor,x,varianceall)})


plot(1:1000,Z[1:1000],type='l', ylab='Z-score', xlab='Genome position')
plot(1000:2000,Z[1000:2000],type='l', ylab='Z-score', xlab='Genome position')
plot(2000:3000,Z[2000:3000],type='l', ylab='Z-score', xlab='Genome position')
plot(3000:4000,Z[3000:4000],type='l', ylab='Z-score', xlab='Genome position')
plot((4000:5000)*10,Z[4000:5000],type='l', ylab='Z-score', xlab='Genome position')
plot(1:length(Z),Z,type='l', ylab='Z-score', xlab='Genome position')



# names(pdfFonts())
plot(control$V1,control$V2,type='l', xlab='Genome position', ylab='Enrichment', ylim=c(0,75))
lines(treat$V1,treat$V2,col=rgb(255,0,0,100,maxColorValue=255))

