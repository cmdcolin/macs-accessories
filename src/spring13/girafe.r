library(ShortRead)
library(girafe)

ex1<-readAligned('S96ChIP_rep1u_eland_result.txt',type="SolexaResult")
ex2<-readAligned('S96ChIP_rep2u_eland_result.txt',type="SolexaResult")

exI1<-as(ex1,"AlignedGenomeRead")






