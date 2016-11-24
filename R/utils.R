#z3_personf depende of this
#pats: matrix with patterns and frcuencie os each pattern
indexPat <- function(data,pats){
  comprimData <- apply(data,MARGIN=1,FUN=paste,collapse="/")
  comprimPats <- apply(pats[,1:ncol(data)],MARGIN=1,FUN=paste,collapse="/")
  index <- lapply(comprimPats,FUN = function(x) {which(x == comprimData)})
  index
}

##the argument are: theta unidmensional and a matrix (z) with the a,d,c parametesr
##calculates the most basic probability in IRT
PR_uni_dico=function(thet,zita){
  
  eta=zita[,1]*thet+zita[,2]
  den=1+exp(-eta)
  
  pr=zita[,3]+((1-zita[,3])/den)
  pr<- ifelse(pr == 0 , sqrt(.Machine$double.eps), pr )  
  pr<- ifelse(pr ==1, 1 - sqrt(.Machine$double.eps), pr) 
  
  
  
}


s_ss <- function(pr,nitems,G){
  sact <- matrix(0,ncol = nitems +1,nrow = G)   
  for(m in 1:G){
    sant <- rep(0,(nitems))
    sant[1] <- 1 - pr[m,1]  
    sant[2] <- pr[m,1]      
    
    for(k in 2:(nitems)){ 
      
      sact[m,1] <- (1-pr[m,k]) * sant[1]   
      
      for(kk in 2:k){ 
        sact[m,kk] <- pr[m,k] * sant[kk-1] + (1 - pr[m,k]) * sant[kk]  
      }
      sact[m,(k+1)] <- pr[m,k] * sant[k] 
      sant <- sact[m,]
    }
    
  }
  return(sact)
}

coll <- function(O, E, mincell = 1){
  for(i in 1L:length(O)){ 
    On <- O[[i]]
    En <- E[[i]]
    
    L <- En < mincell & En != 0
    while(any(L, na.rm = TRUE)){
      if(!is.matrix(L)) break
      whc <- min(which(rowSums(L) > 0L)) 
      if(whc == 1L){ 
        En[2L,] <- En[2L, ] + En[1L,]
        On[2L,] <- On[2L, ] + On[1L,]
        En <- En[-1L,]; On <- On[-1L,]
      } else if(whc == nrow(En)){
        En[nrow(En)-1L,] <- En[nrow(En)-1L, ] + En[nrow(En),]
        On[nrow(On)-1L,] <- On[nrow(On)-1L, ] + On[nrow(On),]
        En <- En[-nrow(En),]; On <- On[-nrow(On),]
      } else {
        ss <- c(sum(On[whc-1L,]), sum(On[whc+1L,]))
        up <- (min(ss) == ss)[1L]
        pick <- if(up) whc-1L else whc+1L
        En[pick,] <- En[pick, ] + En[whc,]
        On[pick,] <- On[pick, ] + On[whc,]
        En <- En[-whc,]; On <- On[-whc,]
      }
      L <- En < mincell & En != 0
    }
    En[En == 0] <- NA
    E[[i]] <- En
    O[[i]] <- On
  }
  return(list("O"=O, "E"=E))
}


x2 <- function(z,patterns,G,FUN){
  nitems <- ncol(patterns$patterns) 
  theta <- patterns$latent_traits 
  frec <- patterns$freq ##
  groups <- quantile(theta,seq(0, 1, length.out = G + 1))
  groups[1] <- groups[1] - 0.1
  groups[G + 1] <- groups[G + 1] + 0.1
  groups.Ind <- findInterval(theta,groups)  
  groups.Ind <- factor(groups.Ind, levels = sort(unique(groups.Ind)))  
  
  thetaG <- tapply(rep(theta, frec), rep(groups.Ind, frec), FUN = FUN) 
  
  
  prs <- matrix(unlist(lapply(thetaG,function(theta){PR_uni_dico(thet=theta,zita = z)})),ncol=nitems,byrow = T)
   # pr=read.table("/home/juan/ppr_ltm.txt")
   # prs=rbind(pr[1:3,],pr[5,],pr[7,],pr[4,],pr[6,],pr[8,],pr[10,],pr[9,])
   # set.seed(17000L)
   # prs=prs+matrix(runif(n = G*nitems,min = 0,max = 0.05),ncol=ncol(prs))
   
  Njs <- as.vector(tapply(frec, groups.Ind, sum))
  Obss2 <- rowsum(frec * patterns$patterns, groups.Ind, reorder = T)/Njs
  
  chi.square <- Njs * (Obss2 - prs)^2/(prs * (1 - prs)) 
  x2 <- colSums(chi.square)
  
  return(x2)
  
}

#################################################################

#prenvelope, envelope_itemf depende of this:
#the arguments are the Specific item parameters and the trait latent 
prenvelope <- function(a,d, cp,  theta){
  D <- 1
  prenve=as.numeric(exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-D*(a*theta+ d)))^(-1))
  prenve=ifelse(prenve<=sqrt(2.2e-16),sqrt(2.2e-16),prenve)
  prenve=ifelse(prenve>=1 - 1e-6,1 - 1e-6,prenve)
  }

#######################################################################
#llikm
#prenvelope, envelope_itemf depende of this:

llikm <- function(params,args){
  a <- params[1]
  d <- params[2]
  c <- params[3]
  
  theta <- args[[1]];
  r <- args[[2]]
  f <- args[[3]]
  item <- args[[4]]
  
  eta <- theta*a+d
  p <- c+(1-c)*(1/(1+exp(-eta)))
  
  p[p<=sqrt(2.2e-16)]=sqrt(2.2e-16)
  p[p >=  1 - 1e-6] =  1 - 1e-6
  
  sum(r[,item]*log(p)+(f[,item]-r[,item])*log(1-p))
  
}

##########################################################
##envelope_itemf depends of thi, is the same function mvtnorm::mvnorm but correcting the eigenvalues

rmvnorm2=function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
          method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE) 
{
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma)) 
    stop("mean and sigma have non-conforming size")
  method <- match.arg(method)
  R <- if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    ev$values=abs(ev$values)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  }
  else if (method == "svd") {
    s. <- svd(sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
  }
  else if (method == "chol") {
    R <- chol(sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% 
    R
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}

######################################################
###normalize: input a matrix, each row is considerted an vector, reurns a list with the norm
###of each vector and the vectors normalized

normalizee=function(matriz){
 
   t(apply(matriz,MARGIN = 1,function(x) {
    norm=sqrt(sum(x*x))
    y=x/norm
    c(norm,y)
    } ))

}



