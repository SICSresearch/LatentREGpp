##################################################################
#
#	CODE FOR CLUSTERING MAIN DIRECTIONS OF UNAL TEST VIA OUR ALGORITHM
#	1) Compute via Multi_Uni + NOHARM values of discrimination
#	2) Normalize
#	3) Cluster via OUR ALGORITHM
#
##################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Necesary Functions (clustering function)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

find_cluster<- function(data, ang=22.5, h7=0.8, q_proy=0.75)
{
  i=1
  aux_data_fix<- list()
  CLUST<- list()
  continue=TRUE
  data_work<- data
  
  
  while(continue==TRUE)
  {
    acp<- dudi.pca(data_work, scannf=FALSE, nf=2)
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
  for(i in 1:length(CLUST_aux))
  {
    if(is.vector(CLUST_aux[[i]])==FALSE){
      CLUST_aux[[j]] = CLUST_aux[[i]]
      acp_clust<- dudi.pca(CLUST_aux[[j]], scannf=FALSE, nf=2)
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