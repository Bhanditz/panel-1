 library(panel)
 sim.ord<-source("simord.dat")[[1]]
 source("simfuns")
 sim.fit<-panel(sim.ord, qfun.ord, gamma.ord, qderiv.ord, 7, 4, 3, T)
