# rm(list=ls())
# setwd("~/Desktop/SICS/simulacion de datasets/")
# source("simulate.poly.multi.R")
# library(FactoMineR)
# 
# simm=simulate.poly.multi(size.cluster = c(25,25,25,25),dim.data = 4,
#                          sample.size = 1000,rep(2,100))
# 
# datos=simm$data
# size.cluster=size.cluster = c(25,25,25,25)

library(FactoMineR)
itemfix=function(datos, size.cluster){

  items=NULL
  
  ind2=cumsum(size.cluster)
  ind1=c(1,ind2[-length(ind2)]+1)
 
   for(i in 1:length(size.cluster)){
     acp=PCA(datos[,ind1[i]:ind2[i]],graph = F)
     items[i]=which(acp$var$cor[,1]==max(acp$var$cor[,1]))
   }
  
return(c(items[1],items[-1]+ind2[-length(ind2)]))
}

