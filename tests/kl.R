# test the data set given in Kalbfleisch and Lawless
 library(panel)
 data(kldata)
 source("klsource")
 klf<-panel(kldata, qfun.kl, theta.kl, qderivs.kl, 8,3,4,T)
