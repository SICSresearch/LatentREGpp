##################################################################
#
#	CODE FOR CLUSTERING MAIN DIRECTIONS OF UNAL TEST VIA OUR ALGORITHM
#	1) Compute via Multi_Uni + NOHARM values of discrimination
#	2) Normalize
#	3) Cluster via OUR ALGORITHM
#
##################################################################
#rm(list=ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("E:/Juan David/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(IRTpp)
library(FactoMineR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Necesary Functions (clustering function)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

find_cluster<- function(data = data7d, ang=22.5, h7=0.8, q_proy=0.75)
{
  i=1
  aux_data_fix<- list()
  CLUST<- list()
  continue=TRUE
  data_work<- data
  
  
  while(continue==TRUE)
  {
    acp<- ade4::dudi.pca(data_work, scannf=FALSE, nf=2)
    names(acp)
    index_aux<- -1*acp$co[,1]
    index_aux2<- acp$co[,2]
    
    which(index_aux== max(index_aux))
    angle<- atan(index_aux2/index_aux)*180/pi
    aux_fix<- which(index_aux== max(index_aux))
    ang_fix<-angle[aux_fix]
    ang_fix_low<- ang_fix - ang
    ang_fix_upp<- ang_fix + ang
    
    
    aux_data_fix[[i]]<- which(angle<ang_fix_upp & angle > ang_fix_low & index_aux >quantile(index_aux, q_proy))
    CLUST[[i]] <- data_work[,(aux_data_fix[[i]])]
    eig_vals<- eigen(cor(data_work))$values
    
    
    data_work<- data_work[, -(aux_data_fix[[i]])]
    if(ncol(data_work)<10) {continue=FALSE}
    if(eig_vals[7]/eig_vals[1]>h7) {continue=FALSE}
    
    i=i+1
    print(i)
  } #END WHILE
  return(clusters=CLUST)
}# END FUNCTION !!!


#**************************************************
#	Principal axes function
#**************************************************

Principal_axes<- function(CLUSTERS, as_matrix=TRUE)
{
  j=1
  CLUST_aux<- CLUSTERS
  P_axes_filter<- list()
  for(i in 1:length(CLUST))
  {
    if(is.vector(CLUST[[i]])==FALSE){
      CLUST_aux[[j]] = CLUST[[i]]
      acp_clust<- ade4::dudi.pca(CLUST_aux[[j]], scannf=FALSE, nf=2)
      P_axes_filter[[j]]<- -acp_clust$li[,1]
      j=j+1
    }
    
  }	
  if(as_matrix==TRUE){P_axes_filter<-do.call(cbind,P_axes_filter)}
  
  return(P_axes_filter)
}# END FUNCTION PRINCIPAL AXES!!!


#AAA<- Principal_axes(CLUST)





#**************************************************
#	Reclassification function
#**************************************************

reclass_data_Praxes<- function(data, Principal_axes)
{
  P_axes_filter<- Principal_axes
  CLUST_REC_FILTER<- list()
  clust_reclass_filter<- vector()
  
  
  for(i in 1:ncol(data))
  {
    proy<- t(P_axes_filter)%*%data[,i]
    clust_reclass_filter[i]<- which(proy == max(proy))[1]
    #print(which(proy == max(proy)))
  }
  
  for(i in 1:length(unique(clust_reclass_filter)))
  {
    
    CLUST_REC_FILTER[[i]]<- data[, clust_reclass_filter==i]
    head(CLUST_REC_FILTER[[i]])
  }
  
  return(CLUST_REC_FILTER)
}# END reclass_data_Praxes function!!!

find_temporal_cluster = function(p,d) {
  pitem <- p
  ddim <- d
  vector = NULL
  for (i in ddim:1) {
    tmp = (pitem+ddim-1)%/%ddim;
    ddim = ddim - 1
    pitem = pitem - tmp
    vector = c(vector,tmp)
  }
  return(vector)
}

#AAA<- Principal_axes(CLUST)
#CLUSTB<- reclass_data_Praxes(data=data, Principal_axes=AAA)
#unlist(lapply(CLUSTB, ncol))



#
#
#	TEST 1: Data with noise
#
#

#data<- read.table("data_without_NA.txt", header=T, sep=";")[,-1]
#head(data); dim(data)
#sc<- c(15,26,29,29,14)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Compute Inivals
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#!!!!! STOP!!!!!
# first charge the Inivals_MultiUni_NOHARM function
# look for "All_FUNCTIONS_multi_uni_plus_NOHARM_IRTpp.R" archive
#!!!!! STOP!!!!!

#III<- inivals_MultiUni_NOHARM(data=data, size.cluster=sc, model="3PL", find.restrictions=TRUE,verbose=TRUE, probit=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Normalize Discrimination vectors per item
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#A_matrix<- III$coefs[,1:5]  ##length size cluster
#Sigma<- cov(A_matrix)
#CholSigm<- chol(Sigma)  
#t(chol(corr))%*%chol(corr)
#t(chol(Sigma))%*%chol(Sigma)
#A_asterisco<- A_matrix%*%solve(CholSigm)
#round(cov(A_asterisco),10)
#Betas<- normalize(t(A_asterisco))
# colSums(Betas^2)
#cov(t(Betas))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#++++++ Find Clusters+++++++#
#CLUST<- find_cluster(data=t(A_matrix), ang=22.5, h7=0.9, q_proy= 0.6)
#unlist(lapply(CLUST,ncol))
#paste("N clust",length(CLUST))


#acpc<-PCA(CLUST[[1]])
#barplot(acpc$eig[,1])
#++++++ Find Reclassify+++++++#

#P_axes_filter<- Principal_axes(CLUST)
#CLUSTB<- reclass_data_Praxes(data=t(A_matrix), Principal_axes=P_axes_filter)
#unlist(lapply(CLUSTB, ncol)); paste("N clust",length(CLUSTB))











#
#
#	TEST 2: Data without noise
#
#
