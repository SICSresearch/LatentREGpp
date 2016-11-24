
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

