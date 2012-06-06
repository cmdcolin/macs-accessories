loaddir<-function(name) {
  oldpath=getwd()
  setwd(name)
  env<-new.env() 
  s96=WiggleClass('S96', env)
  hs959=WiggleClass('HS959', env)
  assign('s96',s96,env) 
  assign('hs959',hs959,env)
  s96$loadWiggles() 
  hs959$loadWiggles()
  ###
  assign('s96bed',read.table('S96/S96_peaks.bed'),env)
  assign('s96overlap',read.table('S96/S96_overlap.bed'),env)
  assign('s96unique',read.table('S96/S96_unique.bed'),env)
  assign('hs959bed',read.table('HS959/HS959_peaks.bed'),env)
  assign('hs959overlap',read.table('S96/S96_overlap.bed'),env)
  assign('hs959unique',read.table('S96/S96_unique.bed'),env)
  # Match indexes 
  id1=match(get('s96overlap',env)$V4,get('s96bed',env)$V4) 
  id2=match(get('s96unique',env)$V4,get('s96bed',env)$V4)
  id3=match(get('hs959overlap',env)$V4,get('hs959bed',env)$V4) 
  id4=match(get('hs959unique',env)$V4,get('hs959bed',env)$V4)
  assign('id1',id1,env)
  assign('id2',id2,env)
  assign('id3',id3,env)
  assign('id4',id4,env)
  assign('name', name)
  setwd(oldpath)
  env
}
e1=loaddir('macs1.4.1')
e2=loaddir('macs1.4.2')
e3=loaddir('macs-bad-run')

plotTotalReads<-function(env, t) {
  s96=get('s96',env)
  hs959=get('hs959',env)
  r1=s96$getTotalReads(get('s96bed',env)) 
  r2=hs959$getTotalReads(get('s96bed',env))
  assign('r1',r1,env)
  assign('r2',r2,env)
  id1=get('id1',env)
  id2=get('id2',env)
  print(tail(sort(r1)))
  print('----')
  print(tail(sort(r2)))
  plot(r1,r2,ylab='HS959 reads',xlab='S96 reads',pch='*') 
  points(r1[id1],r2[id1],pch=1,col='pink') 
  points(r1[id2],r2[id2],pch=1,col='green') 
  title(t)
}
plotTotalReads(e1,'Total S96 peak reads vs HS959 synteny MACS 1.4.1 bw=25')
plotTotalReads(e2,'Total S96 peak reads vs HS959 synteny MACS 1.4.2 bw=25')
plotTotalReads(e3,'Total S96 peak reads vs HS959 synteny MACS 1.4.1 bw=100 [Incorrect]')


debug=TRUE
plotMaxAvgReadsS96<-function(env, t) {
  s96=get('s96',env)
  hs959=get('hs959',env)
  ###########
  # Use Max avg reads over windows
  rma1=s96$getMaxAvgReads(get('s96bed',env),100)
  rma2=hs959$getMaxAvgReads(get('s96bed',env),100)
  id1=get('id1',env)
  id2=get('id2',env)
  assign('rma1',rma1,env)
  assign('rma2',rma2,env)
  ###################### 
  plot(rma1,rma2,ylab='Max Avg HS959 reads',xlab='Max Avg S96 reads',pch='*')
  points(rma1[id1],rma2[id1],pch=1,col='red')
  points(rma1[id2],rma2[id2],pch=1,col='green')
  title('Max Average reads S96 peaks  w=100')
  legend('bottomright', legend=c('shared', 'unique'), fill=c('red', 'green'))
  plot(rma1,rma2,ylab='Max Avg HS959 reads',xlab='Max Avg S96 reads',pch='*',xlim=c(10,30),ylim=c(1,20))
  points(rma1[id1],rma2[id1],pch=1,col='red')
  points(rma1[id2],rma2[id2],pch=1,col='green')
  title('Max Avg reads S96 peaks w=100 (Zoom)')
  legend('bottomright', legend=c('shared', 'unique'), fill=c('red', 'green'))
  
}
plotMaxAvgReadsS96(e2)

plotTotalReadsHS959<-function(env) {
  s96=get('s96',env)
  hs959=get('hs959',env)
  rma3=hs959$getMaxAvgReads(get('hs959bed',env),100) 
  rma4=s96$getMaxAvgReads(get('hs959bed',env),100)
  assign('rma3',rma3,env)
  assign('rma4',rma4,env)
  id3=get('id3',env)
  id4=get('id4',env)
  
  plot(rma3,rma4,ylab='Max Avg S96 reads',xlab='Max Avg HS959 reads',pch='*')
  points(rma3[id3],rma4[id3],pch=1,col='blue')
  points(rma3[id4],rma4[id4],pch=1,col='red')
  title('Max Average reads HS959 peaks w=100')
  legend('bottomright', legend=c('shared', 'unique'), fill=c('blue', 'red'))
  
  # 
  #Zoom 
  plot(rma3,rma4,ylab='Max Avg S96 reads',xlab='Max Avg HS959 reads',pch='*',xlim=c(6,30),ylim=c(1,31))
  points(rma3[id3],rma4[id3],pch=1,col='blue')
  points(rma3[id4],rma4[id4],pch=1,col='red')
  title('Max Average reads HS959 peaks w=100 (Zoom)')
  legend('bottomright', legend=c('shared', 'unique'), fill=c('blue', 'red')) 
}
plotTotalReadsHS959(e2)