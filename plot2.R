s96treat = read.table('S96/S96_MACS_wiggle/treat/S96_treat_afterfiting_chr01.fsa.wig.gz', skip=2)
s96control = read.table('S96/S96_MACS_wiggle/control/S96_control_afterfiting_chr01.fsa.wig.gz', skip=2)
hs959treat= read.table('HS959/HS959_MACS_wiggle/treat/HS959_treat_afterfiting_chr01.fsa.wig.gz', skip=2)
hs959control= read.table('HS959/HS959_MACS_wiggle/control/HS959_control_afterfiting_chr01.fsa.wig.gz')

s96bed=read.table('S96/S96_peaks.bed')
for(i in 1:length(s96bed$V1)) {
  start=s96bed$V2[i];
  end=s96bed$V3[i];
  peak=s96treat[s96treat$V1>start & s96treat$V1<end,];
  plot(1:length(peak$V2), peak$V2,type='l');
  totalreads=apply(peak$V2,1,sum);
  print(totalreads)
}