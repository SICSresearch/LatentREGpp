#
#
# CALCULATE INIVALS MULTI UNI
# 
#
#
#'@name inivals_MultiUni
#'@title Initial Values for one or more dimention models
#'@description Find initial values for any data and model
#'@param data The end is near
#'@param size.cluster waka said
#'@param model and the fools cry
#'@param find.restrictions because why not
#'@param verbose waka waka
#'@export
inivals_MultiUni<- function(data, size.cluster, model="2PL",find.restrictions=FALSE, verbose=FALSE){
  #Input:
  # data: (matrix) dichotomous data set
  # size.cluster: (vector) size of each cluster whom will be adjusted a unidimensional model
  # find.restrictions: (boolean) if find.restrictions, an algorithm decides which vectors are restricted, otherwise first vector of d clusters are restricted 
  # verbose: (boolean) print message in several stages of the processs
  
  #Output:
  # coefs: (matrix) Initial values for the specified model computed via MultiUni models
  # corr: (matrix)  Correlations between constructs 
  
  if(verbose) print("Start process to compute Multidimensional Inivals via Unidimensinal Models.")
  
  #separate cluster
  cluster<- list()
  start = 1
  for(i in 1:length(size.cluster)){
    cluster[[i]]<- data[,start:(start + (size.cluster[i] - 1 ))]
    start = start + size.cluster[i] 
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Reorganize clusters if necesary to find items to
  #restrict
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(find.restrictions)
  {
    
    for(i in 1:length(size.cluster)){
      C<- cor(cluster[[i]])
      E<- RSpectra::eigs(C, 1)   
      Var_coords<- E$vectors*sqrt(E$values)   
      axe.1<- which(Var_coords==max(Var_coords))
      aux_names<-  colnames(cluster[[i]])
      cluster[[i]]<- cbind(cluster[[i]][,axe.1], cluster[[i]][,-axe.1])
      colnames(cluster[[i]])<- c(aux_names[axe.1], aux_names[-axe.1])
    }
  } # END REORGANIZE CLUSTERS
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Recover Unidimensional IRT models to obtain
  #Initial values for NOHARM
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fit_latentregpp<- list()
  traits<- list()
  coef<- list()
  data.traits<- list()
  
  for(i in 1:length(size.cluster)){ # == for(i in 1:nc){
    
    fit_latentregpp[[i]]<- latentreg(data=cluster[[i]], dim=1, model=model, save_time=TRUE, verbose = FALSE)
    z <- fit_latentregpp[[i]]$zetas
    coef[[i]]<- z
    
    traits[[i]]<- ltraits(data= cluster[[i]], dim=1, model="2PL", zetas =fit_latentregpp[[i]]$zetas,  init_traits = NULL, method = "EAP")
    
    
    pattern.matrix<- traits[[i]]$latent_traits[,1]
    pattern.traits<- traits[[i]]$latent_traits[,1]
    data.traits[[i]]<- traits[[i]]$latent_traits[,1]

    
    
    if(verbose) {print(paste("Done", i, "of",length(size.cluster), "dimensions via MultiUni", sep=" " ))}
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Transform "b" (difficultie parameter) into "d"
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  aux1<-lapply(coef,function(x) -x[,2]*x[,1])
  
  #  coef<- lapply(fit_ltm, coef)
  for(i in 1:length(size.cluster)){
    coef[[i]][,2] <- aux1[[i]]
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # CUT EXTREME VALUES     
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for(i in 1:length(size.cluster))
  {
    coef[[i]][,3]<- ifelse(coef[[i]][,3] > 0.3, 0.3, coef[[i]][,3])
    coef[[i]][,2]<- ifelse(coef[[i]][,2] > 3, 3, coef[[i]][,2])
    coef[[i]][,2]<- ifelse(coef[[i]][,2] < -3, -3, coef[[i]][,2])
    coef[[i]][,1]<- ifelse(coef[[i]][,1] > 3, 3, coef[[i]][,1])
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Found Correlation Matrix and Transform A into A^*
  #   and build d and c
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  theta<- matrix(NA, ncol = length(size.cluster), nrow = nrow(data))
  coef_d = coef_c = c()
  A<- list()
  for(i in 1:length(size.cluster)){
    theta[,i]<- data.traits[[i]]
    A[[i]]<- coef[[i]][,1]
    coef_d = c(coef_d, coef[[i]][,2])
    coef_c = c(coef_c, coef[[i]][,3])
  }
  
  aux2<- matrix(NA, ncol=length(size.cluster), nrow= sum(size.cluster))
  for(i in 1:length(size.cluster)){
    col=i
    if(i==1){
      aux2[,col] <- c(as.numeric(c(A[[i]])), rep(NA, length(aux2[,1]) -  length(c(A[[i]]))))
    }else if(i == length(size.cluster)){
      aux2[,col] <- c(rep(NA, length(aux2[,1]) -  length(c(A[[i]]))), as.numeric(c(A[[i]])))
    }else if(i != 1 &&  i!= length(size.cluster)){
      aux2[,col]<- c(rep(NA, sum(size.cluster[1:(col-1)])),
                     as.numeric(c(A[[i]])),
                     rep(NA, sum(size.cluster[(col+1):length(size.cluster)])))
    }
  }
  
  A_matrix<- ifelse(is.na(aux2), 0,aux2)
  sigma<- cov(theta)
  corr<- cor(theta)
  # corr<- matrix(c(1,.5,.5,1), ncol=2)
  # sigma<- matrix(c(.9,.8,.8,.9), ncol=2)
  
  Cholcorr<- chol(corr) 
  CholSigm<- chol(sigma)
  
  #t(chol(corr))%*%chol(corr)
  #t(chol(sigma))%*%chol(sigma)
  
  A_asterisco<- A_matrix%*%CholSigm%*%solve(Cholcorr)
  A_asterisco[is.na(aux2)] = 0
  
  #     A%*%t(theta)
  #     A_asterisco%*%Cholcorr%*%solve(CholSigm)%*%t(theta)
  
  if(verbose) print("Done")
  list(coefs = cbind(A_asterisco, coef_d, coef_c), corr = corr, uni.traits=theta)
  
} #END FUNCTION .... inivals_MultiUni

#
#
# CALCULATE INIVALS MULTI UNI + NOHARM
# (complete Process)
#
#

#'@name inivals_MultiUni_NOHARM
#'@title Initial Values with Noharm for one or more dimention models
#'@description Initial Values with Noharm for one or more dimention models
#'@param data The end is near
#'@param size.cluster waka said
#'@param model and the fools cry
#'@param find.restrictions because why not
#'@param verbose waka waka
#'@param probit This is a waka test, waka
#'@export
inivals_MultiUni_NOHARM<- function(data, size.cluster, model="2PL", find.restrictions=FALSE,  correlated= FALSE,verbose=FALSE, probit=FALSE)
{
  suppressWarnings({
    fit_uni = inivals_MultiUni(data=data, size.cluster=size.cluster, find.restrictions=find.restrictions,verbose=verbose)
  })
  fit = inivals_NOHARM(dat = data, init_uni = fit_uni$coefs, dim_clust = size.cluster, corr = fit_uni$corr, verbose = verbose, probit = probit)
  
  if(correlated==FALSE) 
  {
    A = fit$coefs[,1:length(size.cluster)]
    C = cov(A); E= eigen(C); V=E$vectors; D= diag(E$values); S= V%*%sqrt(D)%*%t(V)
    A_ast<- t(solve(S)%*%t(A))
    colnames(A_ast) = colnames(A)
    fit$coefs[,1:length(size.cluster)]<- A_ast
  }
  fit
}


#III<- inivals_MultiUni_NOHARM(data=DATOS, size.cluster=size.cluster, model="2PL", verbose=TRUE)

#fit_uni = inivals_MultiUni(data=DATOS, size.cluster=size.cluster, find.restrictions=FALSE, verbose=TRUE)
#inivals_NOHARM(dat = data, init_uni = fit_uni$coefs, dim_clust = size.cluster, corr = fit_uni$corr, verbose = TRUE, probit = FALSE)











#
#
# CALCULATE INIVALS NOHARM
# 
#
#




#'@name inivals_NOHARM
#'@title Noharm's initial values
#'@description Auxiliar function, maybe
#'@param dat waka
#'@param init_uni waka waka
#'@param dim_clust waka wakiri waka waka
#'@param corr wukuru waka wa
#'@param verbose kawa waka burm bum
#'@param probit waka false waka
#'@export
inivals_NOHARM<- function(dat, init_uni, dim_clust, corr, verbose = FALSE, probit = FALSE){
  #Input:
  # dat: (matrix) it's dichotomic dataset
  # init_uni: (matrix) it's initial values for coefs. 
  # dim_clust: array with clusters dimensions
  # corr: (matrix) it's intial value of correlation matrixmatrix
  # Verbose: (boolean) print message in start and finish process
  # probit: (boolean) if probit = False then the model is logit else the model is probit
  
  #Output:
  # coefs: (matrix) initial values to multidimensional IRT model
  # corr: (matrix) correlations between constructs 
  
  
  if(verbose) print("Start process to calculate multidimensional init vals.")
  
  #Control section
  if(length(unique(dim(corr))) > 1 || !isSymmetric(corr) || sum(corr > 1) > 0) stop("The matrix corr is not correlation matrix")
  if(!is.logical(probit)) stop("The value of probit is not logical")
  if(sum(as.matrix(dat) %in% c(0,1)) + sum(is.na(dat)) != dim(dat)[1] * dim(dat)[2]) stop("The matrix dat is not binary")
  if(is.null(dim_clust) | length(dim_clust) != dim(corr)[1]) stop("The value of dim_clust is invalid")
  if(dim(init_uni)[2] - 2 != dim(corr)[2]) stop("Dimensions of init_uni matrix and corr matrix do not correspond")
  
  #Convert to probit
  if(!probit){
    sub_init = init_uni[,c(1:(ncol(init_uni) - 1))]
    sub_init = sub_init / 1.702
    init_uni[,c(1:(ncol(init_uni) - 1))] = sub_init
  }
  
  #Create Fpatt and Fval
  dim = dim(corr)[2]
  Fpatt = matrix(1, nrow = ncol(dat), ncol = dim)
  start_act = 1
  for(i in 1:length(dim_clust)){
    Fpatt[start_act,] = rep(0, dim)
    start_act = start_act + dim_clust[i]
  }
  
  Fval = init_uni[,1:dim]
  
  
  
  #Created Ppatt
  Ppatt = matrix(0, ncol = dim, nrow = dim)
  diag(Ppatt) = 0
  
  #Get noharm fit
  fit = sirt::noharm.sirt( dat = dat, Ppatt = Ppatt, Fpatt = Fpatt, Fval = Fval, Pval = corr)
  
  #build list for return
  if(!probit){
    coef_a = fit$loadings * 1.702
    coef_b = fit$final.constants * 1.702
  }else{
    coef_a = fit$loadings
    coef_b = fit$final.constants
  }
  
  coefs = cbind(coef_a, coef_b, init_uni[,ncol(init_uni)])
  colnames(coefs) = c(paste("a_", 1:dim, sep = ""), "d", "c")
  corr = fit$factor.cor
  colnames(corr) = c(paste("a_", 1:dim, sep = ""))
  rownames(corr) = c(paste("a_", 1:dim, sep = ""))
  ret = list(coefs = coefs, corr = corr)
  
  if(verbose) print("Done")
  
  ret
} #END FUNCTION .... inivals_NOHARM



#inivals_MultiUni_NOHARM(data=datam, size.cluster=size.cluster, find.restrictions=FALSE, verbose=TRUE, probit=FALSE)





####################################################################
#
# FUNCTION TO FIND LATENT TRAITS FOR ALL INDIVIDUALS
# HAVING LATENT TRAIT FOR PATTERNS
#
####################################################################
#'@name traits_patt2data
#'@title Returns a list with the trait of every individual
#'@description Returns a list with the trait of every individual
#'@param pattern Is the list with all response patterns
#'@param pattern.traits Is the list with the trait for each pattern response
#'@param the data used to make the estimation of traits and paramters
#'@export
traits_patt2data<- function(pattern, pattern.traits, data)
{
  PATTERNS <- pattern
  DATA<- data
  
  D<-as.matrix(dist(rbind(as.matrix(PATTERNS), as.matrix(DATA))))
  D_real<-D[1:nrow(PATTERNS),(nrow(PATTERNS)+1):nrow(D)]
  
  P_return<- list()
  for(i in 1:ncol(D_real))
  {
    P_return[[i]]<-as.numeric(which(D_real[,i]==0))
  }
  
  data.traits<- vector()
  for(i in 1:length(P_return))
  {
    data.traits[i]<- pattern.traits[P_return[[i]]]
  }
  
  return(data.traits)
}


#
#
#
traits_by_individuals<- function(data, patterns)
{
  all_patterns = new.env()
  for ( i in 1:nrow(data) ) {
    key = paste(data[i,], collapse = "")
    all_patterns[[key]] = append(all_patterns[[key]], i)
  }
  
  traits_by_individuals = list()
  
  for ( l in 1:nrow(patterns) ) {
    key = paste(patterns[l,-ncol(patterns)], collapse = "")
    indexes_pattern_l = all_patterns[[key]]
    for ( h in 1:length(indexes_pattern_l)) {
      index = indexes_pattern_l[[h]]
      traits_by_individuals[[index]] = patterns[l,ncol(patterns)]
    }
  }
  all_individual_traits<-unlist(traits_by_individuals)
  return(all_individual_traits)
}

normalize<- function(y)
{
  m.y<- apply(y,2,function(x) sqrt(sum(x^2)))
  y.normalize<- t(t(y)/m.y)
  return(y.normalize)
}