library(DiffBind)
K27ac_B<-dba(sampleSheet = "K27ac_B.csv",dir=system.file("extra", package="DiffBind"))
K27ac_B_consensus<-dba.peakset(K27ac_B,bRetrieve = TRUE)
write.csv(K27ac_B_consensus,file="C:/Users/Jerry Liu/Documents/K27ac_B_consensus.csv")