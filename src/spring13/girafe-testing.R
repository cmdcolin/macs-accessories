


library(ShortRead)
ex<-readAligned('S96ChIP_rep1u_eland_result.txt',type='SolexaResult')
ex2<-readAligned('S96ChIP_rep2u_eland_result.txt',type='SolexaResult')


exAI <- as(ex, "AlignedGenomeIntervals")
exBI <- as(ex2, "AlignedGenomeIntervals")
table(seq_name(exAI))
table(seq_name(exBI))