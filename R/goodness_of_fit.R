#######################################################################
#' @name z3_personf
#' @title Z3 Person fit statistic
#' @description Calculates the values of statistical Z3 for individuals.
#' @param data a data frame or a matrix with the test.
#' @param zita Matrix whose columns are the estimations of the parameters of the items (discrimination,d=-discrimnation*difficulty, guessing).
#' @param patterns list of patterns response (patterns$patterns), the frequency of each pattern (patterns$freq) and the latent traits (patterns$latent_traits).
#' @export
#'
#' @seealso
#' \code{\link{z3_itemf}}, \code{\link{orlando_itemf}}
#'
#' @references
#'
#' Fritz Drasgow, Michael V. Levine and Esther A. Williams (1985). Appropiateness measurement with polychotomous item response models and standarized indices.
#'
z3_personf <- function(data,zita,patterns){
  
  scores=cbind(patterns$patterns,patterns$latent_traits)
  ninds <- nrow(data) 
  nitems=ncol(data)
  index <- indexPat(data = data,pats = cbind(patterns$patterns,patterns$freq)) 
  scoresTot <- numeric(nrow(data)) 
  for(mm in 1:nrow(scores)){
    scoresTot[index[[mm]]] <- scores[mm,ncol(data) +1]  
  }
  
  P <- lapply(scoresTot,FUN=function(x){PR_uni_dico(thet=x,zita=zita)}) 
  P <- matrix(unlist(P),ncol=nitems,byrow=T)
  dim(P)
  
  LL <- matrix(0,ncol = ncol(P),nrow = nrow(P)) 
  LL[data == 1] <- P[data == 1] 
  LL[data == 0] <- 1 - P[data == 0] 
  LL <- rowSums(log(LL)) 
  
  mu <- sigmaCuad <- rep(0,ninds)
  for( i in 1:nitems){
    Pi <- cbind(P[,i],1 - P[,i])
    logPi <- log(Pi)
    mu <- mu+rowSums(Pi * logPi)
    sigmaCuad <- sigmaCuad + (Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))
  }
  
  Z3 <- (LL - mu) / sqrt(sigmaCuad)
  Z3
  
}

#######################################################################
#' @name z3_itemf
#' @title Z3item fit statistic
#' @description Calculates the values of statistical Z3 for items.
#' @param data a data frame or a matrix with the test.
#' @param zita Matrix which columns are the estimations of the parameters of the items (discrimination,d=-discrimnation*difficulty, guessing).
#' @param patterns list of patterns response (patterns$patterns), the frequency of each pattern (patterns$freq) and the latent traits (patterns$latent_traits).
#' @export
#'
#' @seealso
#' \code{\link{z3_personf}}, \code{\link{orlando_itemf}}
#'
#' @references
#'
#' Fritz Drasgow, Michael V. Levine and Esther A. Williams (1985). Appropiateness measurement with polychotomous item response models and standarized indices.
#'
z3_itemf <- function(data,zita,patterns){
  
  scores=cbind(patterns$patterns,patterns$latent_traits)
  ninds <- nrow(data) 
  nitems=ncol(data)
  index <- indexPat(data = data,pats = cbind(patterns$patterns,patterns$freq)) 
  scoresTot <- numeric(nrow(data)) 
  for(mm in 1:nrow(scores)){
    scoresTot[index[[mm]]] <- scores[mm,ncol(data) +1]  
  }
  
  P <- lapply(scoresTot,FUN=function(x){PR_uni_dico(thet=x,zita=zita)}) 
  P <- matrix(unlist(P),ncol=nitems,byrow=T)
  dim(P)
  
  LL <- matrix(0,ncol = ncol(P),nrow = nrow(P)) 
  LL[data == 1] <- P[data == 1] 
  LL[data == 0] <- 1 - P[data == 0] 
  LL <- colSums(log(LL)) 
  
  mu <- sigmaCuad <- rep(0,nitems)
  for( i in 1:nitems){
    Pi <- cbind(P[,i],1 - P[,i])
    logPi <- log(Pi)
    mu[i] <- sum(Pi * logPi)
    sigmaCuad[i] <- sum(Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))
  }
  
  Z3 <- (LL - mu) / sqrt(sigmaCuad)
  Z3
  
}

#######################################################################
#' @name orlando_itemf
#' @title Orlando's statistic
#' @description Calculate the values of the statistics S_x2 from
#' Maria Orlando and David Thisen (2000).
#' @param patterns list of patterns response (patterns$patterns), the frequency of each pattern (patterns$freq) and the latent traits (patterns$latent_traits).
#' @param G number of quadrature points.
#' @param zeta matrix of estimations of the parameters of the items (alphas, d's, guessing).
#' @param model type of model ( "1PL", 2PL", "3PL" ).
#' @return Orlando's statistic, degrees of freedom and p-value for each item.
#'
#' @seealso
#' \code{\link{z3_itemf}}
#' 
#' @references
#' Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for dichotomous item
#' response theory models. \emph{Applied Psychological Measurement, 24}, 50-64.
#' @export
orlando_itemf <- function(patterns,G,zeta,model){
  zita = zeta;
  
  if(model=="3PL"){mo=3}
  if(model=="2PL"){mo=2}
  if(model=="1PL"){mo=1}
  
  #########
  
  pats <- as.matrix(cbind(patterns$patterns,patterns$freq))
  frec <- patterns$freq  
  patsSinFrec <- patterns$patterns
  nitems <- ncol(patsSinFrec)
  
  seq <- seq(-6,6,length=G)
  pesos <- dnorm(seq)/sum(dnorm(seq))
  Cuad <- matrix(c(seq,pesos),byrow=F,ncol=2)
  
  theta <- Cuad
  w.cuad <- theta[,2] 
  thetaG <- theta[,1] 
  
  pr <- lapply(thetaG,FUN=function(x){PR_uni_dico(thet=x,zita=zita)}) #probabilidad para cada punto de cuadratura
  pr <- matrix(unlist(pr),ncol=nitems,byrow=T)
  
  score <- rowSums(patsSinFrec )  
  Nk <- NULL
  for(i in 1:nitems - 1){
    inds <- which(score == i) 
    patsCoin <- pats[inds,]   
    if(class(patsCoin) == "matrix"){
      if(dim(patsCoin)[1] != 0){     
        Nk[i] <- sum(patsCoin[,ncol(patsCoin)])
      }else{
        Nk[i] <- 0
      }
    }else{
      if(class(patsCoin) == "numeric"){
        Nk[i] <- patsCoin[length(patsCoin)]
      }
    }
  }
  
  O <- list()
  Oik <- matrix(0,ncol = nitems -1 ,nrow = nitems)  
  for(i in 1:nitems - 1){  
    inds <-which(score == i); 
    for(j in 1:(nitems)){  
      patsCoin <- pats[inds,]   
      if(class(patsCoin) == "matrix"){
        if(dim(patsCoin)[1] != 0){     
          Oik[j,i] <- sum(apply(X = patsCoin,MARGIN = 1,FUN = function(x){ifelse(x[j] == 1,yes = x[nitems + 1],0)}))
        }else{
          Oik[j,i] <- 0
        }
      }else{
        if(class(patsCoin) == "numeric"){
          Oik[j,i] <- ifelse(patsCoin[j] == 1,yes = patsCoin[nitems + 1],0)
        }
      }
      O[[j]] <- cbind(Nk-Oik[j,],Oik[j,])
    }
  }
  
  sact <- s_ss(pr,nitems=nitems,G=G)
  Denom <- colSums(matrix(rep(w.cuad,nitems -1 ),ncol = nitems - 1) * sact[,-c(1,ncol(sact))])
  
  nitems <- nitems - 1
  smono <- list()
  
  for(p in 1:(nitems+1)){
    smono[[p]] <- s_ss(pr[,-p],nitems=nitems,G=G)
  }
  
  nitems <- ncol(patsSinFrec)
  E <- list()
  for(i in 1:length(smono)){
    E[[i]] <- cbind(1-colSums(smono[[i]][,-ncol(smono[[i]])]*(pr[,i]*w.cuad))/Denom,colSums(smono[[i]][,-ncol(smono[[i]])]*(pr[,i]*w.cuad))/Denom)
  }
  
  S_X2 <- NULL
  df.S_X2 <- NULL
  for (i in 1:nitems) E[[i]] <- E[[i]] * Nk
  coll <- coll(O, E, mincell = 1)
  O <- coll$O
  E <- coll$E
  for (i in 1:length(O)) {
    S_X2[i] <- sum((O[[i]] - E[[i]])^2/E[[i]], na.rm = TRUE)
    df.S_X2[i] <- sum(!is.na(E[[i]])) - nrow(E[[i]]) - mo
  }
  
  pval <- pchisq(S_X2,df.S_X2,lower.tail = F)
  
  lista <- cbind("S_X2"=S_X2,"df.SX2"=df.S_X2,"p.val.S_X2"=pval)
  
  return(lista)
  
}

#######################################################################
#' @name x2_itemf
#' @title Statistical x2.
#' @description Calculates the statistical x2.
#' @usage x2_itemf(zetas, patterns, G, FUN)
#' @param zetas matrix of estimations of the parameters of the items (alphas, d's, guessing).
#' @param patterns list with Patterns, frequencies and traits.
#' @param G the number of groups, by default 10
#' @param FUN It is the function with which the expected probability, by default median 
#' is calculated in each group.
#' @export
x2_itemf=function(zetas,patterns,G = 10,FUN = median){
  x2(z=zetas,patterns=patterns,G=G,FUN=FUN)
}


#######################################################################
#' @name envelope_itemf
#' @title Confidence Envelops for an item
#' @description Graphs confidence bands of an item, to evaluate 
#' the goodness of fit of the model.
#' @param item a number indicating the item to be evaluated.
#' @param numboot number of iterations bootstrap, used to plot the envelopes. By default 100.
#' @param alpha level of significance to plot the envelopes. By default 0.05.
#' @param item.fit object LatentREGpp::itemfit() type.
#' @param data a dataframe or a matrix with the test data.
#' @param seed the seed to fix the random sample. By default 500L.
#' @return plot with the envelopes and the caracteristic curve of the item.
#'
#' @seealso
#' \code{\link{orlando_itemf}}, \code{\link{z3_itemf}}
#'
#' @references
#'
#' David Thissen, Howard Wainer  D. (1990). Confidence Envelopes for Item Response Theory. \emph{Journal of Educational Statistics, Vol 15, No 2}, 113-128.
#' @export
envelope_itemf=function(item, numboot=100, alpha=0.05, item.fit, data,seed=5000L){
  
  r <- item.fit$r
  f <- matrix(rep(item.fit$f,ncol(data)),ncol=ncol(data),byrow=F)
  theta <- item.fit$quad$theta
  params <- item.fit$zetas
  
  nitems <- ncol(r) 
  x <- seq(from = -6,to = 6,by=0.1)
  y <- as.vector(sapply(X = x,FUN = prenvelope,a = params[item,1],d = params[item,2],cp = qlogis(params[item,3]))) 
  ####hasta aqui no hay problemas de restricciones.
  
  inf <- solve(-hessian(func=llikm,x=as.vector(params[item,]),args=list(theta,r,f,item)))
  
  media <- params[item,] 
  set.seed(seed)
  boot <- rmvnorm2(numboot,mean = media,sigma = inf)
  boot[,3] <- ifelse(boot[,3] >= 1,1 - 1e-6,boot[,3])
  boot[,3] <- ifelse(boot[,3] <= 0,sqrt(2.2e-16),boot[,3])
  boot[,1] <- ifelse(boot[,1] <= 0,sqrt(2.2e-16),boot[,1])
  
  envelop <- matrix(0,nrow = length(x),ncol = 2)
  for(i in 1:length(x)){
    trace <- apply(X = boot, MARGIN = 1,
                   FUN = function(puntoBoot){
                     prenvelope(a = puntoBoot[1], d = puntoBoot[2],cp = qlogis(puntoBoot[3]),theta = x[i])
                   })
    trace <- sort(trace) 
    lower <- trace[floor(length(trace) * alpha / 2)] 
    upper <- trace[ceiling(length(trace) * (1 - (alpha / 2)))] 
    envelop[i,] = c(lower,upper)
    
  }
  
  plot(x,y,xlim=c(-6,6),ylim=c(0,1),type="l",main=paste("Item Characteristic Curve",sep=" "),
       xlab = expression(theta),ylab = expression(P(theta)))
  lines(x,envelop[,1],col="red",lty=2)
  lines(x,envelop[,2],col="red",lty=2)
  
}

