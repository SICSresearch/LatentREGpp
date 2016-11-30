#######################################################################
#' @name findim_hcpc
#' @title Discovering Dimension
#' @description Finds the test dimension for a dichotomous logistic model by using Loadings and Factorial Analysis, according to the Paez & Montenegro methodology.
#' @usage findim_hcpc(data, verbose, probit)
#' @param data A binary matrix or dataframe that holds the response data with N individuals (columns) and P items (rows).
#' @param verbose a boolean, if TRUE all procedures are descibed in console. The adjusted eigen values are allways shown.
#' @param probit a boolean, if FALSE correction of probit to logit model is made (multiplying by 1.702 de discrimination vectors).
#' @details Implementation of a tecnique to evaluate the number of latent constructs presented by data.
#' @return L_data_Clust: list of matrices formed by clusters, each matrix has an associated set of binary vectors of items.
#' @return M_data: matrix of the join of objects in L_data_Clust.
#' @return dim_new: numeric, Dimension found by the algorithm.
#' @return dim_old: numeric, Initial Dimension of the algorithm (set by parallel analysis).
#' @return sc: vector, The size of each cluster of L_data_Clust.
#' @examples
#' \dontrun{
#' file = paste(system.file(package = "LatentREGpp"),"/dataset/TUN.txt",sep = "")
#' TUN = data.matrix(read.table(file))
#' findim_hcpc(TUN)
#' findim_hcpc(TUN, verbose = T)
#' }
#' @references Paez S. Montenegro A. and Pardo C. (2017). Principles and Methodology to Find Dimension on Latent Regression Models. Novel based approach.
#' British Journal of Mathematical and Statistical Psychology. (Submitted)
#' @references John L. Horn (1965). A rationale and test for the number of factors 
#' in factor analysis. Psychometrika, Volume 30, Number 2, Page 179.
#' @references Lebart L, Morineau A, Piron M (1997). Statistique Exploratoire Multidimensionnelle. Dunod.
#' @references Reckase M (2009). Multidimensional item response theory. Springer.
#' @seealso
#' \code{\link{itemfit}}
#' @export
findim_hcpc<- function(data, verbose=FALSE, probit = F)
{
  stopifnot(!is.null(data))

  if(is.null(colnames(data)))
  {	
    colnames(data)<- paste("V",1:ncol(data), sep="")
  }
  
  ## a conservative analysis with different result!
  ##############################################
  #	Using Paralell Analysis to 
  # 	found a initial dimension
  ##############################################
  
  if(verbose) {cat("Using Paralel Analysis \n")}
  paran<- try(paran::paran(data, iterations=500, centile=95, graph=T), silent=T)
  if(class(paran) == "try-error")
  {
   cat("\n\n\n")
  stop("There are not eigen values greater than 1, then it is not even 1 dimension")
  }
  
  #*****	Adjusted eigen values		****#
  #*****	Greater than 1 are selected	****#

  dim<- paran$Retained
  acp<-   ade4::dudi.pca(data, nf=dim, scannf = FALSE)
  Axis<- -acp$li
  Prod<- t(Axis)%*%as.matrix(data)
  
  #*****	Wait user to continue		****#
  #*****						****#
  
 	  
  cat("\nPress [enter] to continue: ")
  Sys.sleep(0.01)
  readline()
  
  ##############################################
  #	Selecting items to be fixed
  # 	items with the max projection 
  #	over each principal direction
  #	of the PCA are fixed
  ##############################################
  
  
  index_prod<- vector()
  for(i in 1:dim)
  {
    j=1
    end=FALSE
    index_prod[i]<- which(Prod[i,] == 	apply(Prod,1, sort, decreasing=T)[,i][j])
    if(i > 1)
    {
      while(end==FALSE) {
        l= length(unique(index_prod))
        if(l < i)
        {
          j=j+1
          index_prod[i]<- which(Prod[i,] == 	apply(Prod,1, sort, decreasing=T)[,i][j])
        }
        
        if(l == i)
        {
          end=TRUE
        }
      }
    }
  }
  
  
  
  ##############################################
  #	Preparing to use noharm
  ##############################################
  
  
  
  PPatt<- matrix(0, ncol=dim, nrow=dim)
  FPatt<-  matrix(1, ncol=dim, nrow=ncol(data))
  FPatt[index_prod,]<- matrix(0, ncol=dim, nrow=dim)
  
  
  PVal<- diag(dim)
  FVal<-  matrix(1, ncol=dim, nrow=ncol(data))
  FVal[index_prod,]<- diag(dim)
  
  
  if(verbose){ cat("\n++++  Computing items' directions  ++++")}
  fit <- sirt::noharm.sirt(dat= data, Ppatt= PPatt, Fpatt= FPatt, Pval=PVal, Fval= FVal)
  
  
  
  #*****	Normalizing item slopes		****#
  #*****	to generate directions		****#
  
  
  A<- fit$loadings.theta
  if(!probit){A<- fit$loadings.theta*1.702}
  norm<- sqrt(apply(A^2, 1, sum))
  A_norm<- A/matrix(norm, ncol=dim, nrow=nrow(A))
  
  
  ##############################################
  #	Using HCPC from FactoMineR
  #	to generate groups of items
  ##############################################
  
  
  if(verbose){ cat("\n++++  Making new clusters which generate new dimension  ++++")}
  if(dim == 1)
  {
  cat("\n\n")
  stop("There is just one dimension")
  }

  acp_dir<- FactoMineR::PCA(A_norm, graph=FALSE)
  hcpc_dir<- FactoMineR::HCPC(acp_dir)
  clust<- hcpc_dir$data.clust[, (dim+1)]
  
  #*****	Number of clusters is 		****#
  #*****	The new dimension 		****#
  
  
  dim.new<- length(unique(clust))  
  
  
  #*****	Saving new clusters		****#
  
  
  CLUST<- list()
  for(i in 1:dim.new)
  {
    CLUST[[i]]<- data[,which(clust==i)]
    colnames(CLUST[[i]])<- colnames(data)[which(clust==i)]
  }
  
  if(verbose){ cat("\n++++  Done 'Dim_found_hcpc' ++++")}
  
  res<- list("L_data_Clust" = CLUST, "M_data" = do.call(cbind, CLUST),"dim_new" = dim.new, "dim_old" = dim, 
             "sc" = unlist(lapply(CLUST, ncol)))
  return(res)
}

