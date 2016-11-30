loglik_M=function(betas.item,r,pt.cuad,item,ind1,ind2,ncatg){

  p=prob_gpcm_M(betas.item,pt.cuad,item,ind1,ind2,ncatg)
  r.item=r[,c(ind1[item]:ind2[item])]
  return(-sum(r.item*log(p)))
  
}

prob=function(betas,pt.cuad,nitems,ncatg){

  probs=list()
  for(g in 1:length(pt.cuad)){###para cada nodo:
  
      node=pt.cuad[g]
      p=prob_gpcm(betas = betas,node = node,nitems = nitems,ncatg)
    
  	  p<-lapply(p, function(x) ifelse(x <= sqrt(2.2e-16) , sqrt(2.2e-16), x )  ) 
  	  p<- lapply(p, function(x) ifelse(x >=  1 - 1e-6,  1 - 1e-6 , x) )
  	  probs[[g]]=p
  }
  return(probs)
}

posteriori=function(pt.cuad,w.cuad,betas,pats,ncatg,npats,dat.dic){
	nitems=length(betas)
	Pi=matrix(NA,nrow=length(pt.cuad),ncol=npats)
	p=prob(betas = betas,pt.cuad = pt.cuad,nitems=nitems,ncatg)

	for(l in 1:nrow(pats)){

    dat.dic.pat=dat.dic[[l]]
  
	    for(g in 1:length(pt.cuad)){
	 
	      	p.node=p[[g]]
	      
	      	expon=mapply(function(x,y) x^y,p.node,dat.dic.pat,SIMPLIFY = F)
	      	prod.catgs=lapply(expon,prod)
	      	Pi[g,l]=do.call(prod,prod.catgs)*(w.cuad[g])
	    }
    }
    post=Pi/matrix(colSums(Pi),nrow=nrow(Pi),ncol=ncol(Pi),byrow=T)
    return(post)
}

rr=function(pt.cuad,w.cuad,betas,pats,ncatg,freqs,npats,dat.dic,pats.cod){
	pii=posteriori(pt.cuad = pt.cuad,w.cuad = w.cuad,betas = betas,pats = pats,ncatg = ncatg,npats,dat.dic)
	freq.matrix=matrix(freqs,nrow=length(pt.cuad),ncol=length(freqs),byrow=T)
	(pii*freq.matrix)%*%pats.cod
}

patrones=function(datos){
  pats <- apply(datos, 1, paste, collapse = "/")
  freqs <- table(pats) 
  nfreqs <- length(freqs)
  X <- unlist(strsplit(cbind(names(freqs)), "/"))
  X <- matrix(as.numeric(X), nfreqs, ncol(datos), TRUE)
  return(list("data"=X,"freqs"=as.vector(freqs)))
}

expand <- function(tab,ncatg) {
  Y <- matrix(0, nrow(tab), ncatg*ncol(tab))
  for(i in 1:nrow(tab)){
    moda <- tab[i,]
    n <- length(moda)
    x <- matrix(0, n, ncatg)
    x[(1:n) + n * (moda - 1)] <- 1
    x <- t(x)
    x <- as.vector(x)
    Y[i,] <- x 
  }
  return(Y)
} 

expand.list <- function(pats,ncatg,npats) {
  dat.dicotom=list()
  for(s in 1:npats){
  
	  P<- length(ncatg)
	  Y <- list()
	  
	  for(item in 1:P)
	  {	
		  y <- vector(length=(ncatg[item]))
		  y[pats[s,item]] =T	
		  y<- ifelse(y,1,0)
		  Y[[item]]<- y
	  }
	  dat.dicotom[[s]]=Y
  }
  return(dat.dicotom)
}

prob_gpcm_M=function(betas.item,pt.cuad,item,ind1,ind2,ncatg){
  
  pr=matrix(NA,nrow=length(pt.cuad),ncol=length(c(ind1[item]:ind2[item])))
  
  for(i in 1:nrow(pr)){ 
	  node=pt.cuad[i] 
	  
	  etas=betas.item[length(betas.item)]*(node-betas.item[1:(length(betas.item)-1)])
	  
	  cumsumetas=cumsum(etas)
	  den=1+sum(exp(cumsumetas))
	  pr[i,1]=1/den
	  
	  for(k in 2:(ncatg[item])){ 
	      pr[i,k]=exp(etas[k-1])*pr[i,k-1]
	  }
  }
  pr[pr<=sqrt(2.2e-16)]=sqrt(2.2e-16)
  pr[pr>=1 - 1e-6]=1 - 1e-6
  return(pr)
  
}

prob_gpcm=function(node,betas,nitems,ncatg){

	etas=lapply(betas,function(x,node){x[length(x)]*(node-x[1:(length(x)-1)])},node=node)

	cumsumetas=lapply(etas,function(x)cumsum(x))
	den=lapply(cumsumetas,function(x){1+sum(exp(x))})
	p=list()
	for(i in 1:nitems){
		p[[i]]=vector(length = ncatg[i]);
		p[[i]][1]=1/den[[i]]
	}

	for(i in 1:nitems){ 
	  for(k in 2:(ncatg[i])){ 
	    p[[i]][k]=exp(etas[[i]][k-1])*p[[i]][k-1]
	  }
	}
	p<-lapply(p, function(x) ifelse(x <= sqrt(2.2e-16) , sqrt(2.2e-16), x )  ) 
	p<- lapply(p, function(x) ifelse(x >=  1 - 1e-6,  1 - 1e-6 , x) )

	# p=lapply(p,function(x){x[-1]})
	return(p)
}

#'@name estim_gpc
#'@title Estimation gpc model
#'@description Estimates the test parameters according to the gerenalized parcial credit model (gpc) model.
#'This model is unidimensional only.
#'@param data The matrix containing the answers of tested individuals
#'@examples
#'\dontrun{
#' data = simulate_polytomous$data
#' estim=estim_gpc(data = datos)
#'}
#'@export
estim_gpc=function(data){      
  
  ###numero de categorias por item y el item a estimar.
  ncatg <- apply(data, 2, function (x) if (any(is.na(x))) length(unique(x)) - 1 else length(unique(x)))
  ind1 <- c(1, cumsum(ncatg[-ncol(data)]) + 1)
  ind2 <- cumsum(ncatg)
  nitems=ncol(data)
  
  ###valores iniciales:
  betas=list()
  # param=list()
  for(i in 1:nitems){
     betas[[i]]=c(seq(-1.5,1.5,length=ncatg[i]-1),1)
    #  betas[[i]]=as.vector(valini[i,])
    # # param[[i]]=c(betas[[i]][1],log(diff(betas[[i]][-ncatg[i]])),1)
  }
  
  
  ###nodos y pesos:
  Cuad = gauss.quad(n=40,"hermite") #nodos y pesos
  pt.cuad = Cuad[[1]]*sqrt(2) 
  w.cuad = Cuad[[2]]/sqrt(pi) #transforman los pesos
  
  #patrones y frecuencias
  pats=patrones(data)$data
  npats=nrow(pats)
  freqs=patrones(data)$freqs
  
  #datos expandidos por item
  pats.cod=expand(pats,ncatg)
  
  #datos dicotomizados para cada item, y para todas las categorias
  dat.dic=expand.list(pats=pats,ncatg = ncatg,npats)
  
  
  ########### ALGORITMO EM ###########
  
  #em:
  seguir = TRUE
  mm = 0
  betas.ant=betas
  #param=param.ant=param.true  #deltas verdaderos
  #betas=betas.ant=param.sim  #betas verdaderos
  # warn=NULL
  
  while(seguir) {
    
    ##contando los ciclos:
    mm = mm+1
    
    #PASO E
    r=rr(pt.cuad = pt.cuad,w.cuad = w.cuad,betas = betas,pats = pats,ncatg = ncatg,freqs = freqs,npats=npats,dat.dic=dat.dic,pats.cod=pats.cod)
    print(sum(r))
    
    #PASO M
    for(item in 1:nitems){  
      betas.item=betas[[item]]
      opt<-optimx(par = betas.item,fn = loglik_M,method = "BFGS",r=r,pt.cuad=pt.cuad,item=item,ind1=ind1,ind2=ind2,ncatg=ncatg)
      betas.item=unlist(opt[1:length(betas.item)])
      betas[[item]]=betas.item
    }
    
   # print(mm)
    
    ###evaluando convergencia:
    if(max(mapply(function(x,y){max(x-y)},betas.ant,betas)) < 10^(-3)){
      seguir = FALSE
    }
    betas.ant=betas
    # if(mm==500){warn="no converge";seguir=FALSE}
  }
  retorno=list(pt.cuad=pt.cuad,w.cuad = w.cuad,betas = betas,
               pats = pats,ncatg = ncatg,dat.dic = dat.dic)
  return(retorno)
}

#'@name eap_gpc
#'@title Latent Trait Estimation with EAP for gerenalized parcial credit model (gpc) model 
#'@description Estimates latent traits for gerenalized parcial credit model (gpc) model.
#'Just for unidimensional models.
#'@param estim Output object from the estim_gpc function.
#'@examples
#'\dontrun{
#' data = simulate_polytomous()$data
#' estim=estim_gpc(datos = datos)
#' eap_gpc = eap_gpc(estim)
#' plot(density(eap_gpc))
#'}
#'@export
eap_gpc=function(estim){
  betas=estim$betas
  pt.cuad=estim$pt.cuad
  w.cuad=estim$w.cuad
  pats=estim$pats
  ncatg=estim$ncatg
  dat.dic=estim$dat.dic
  nitems=length(betas)
  npats=nrow(pats)
  denom=matrix(NA,nrow=length(pt.cuad),ncol=npats)
  num=matrix(NA,nrow=length(pt.cuad),ncol=npats)
  p=prob(betas = betas,pt.cuad = pt.cuad,nitems=nitems,ncatg)
  
  for(l in 1:nrow(pats)){  
    
    dat.dic.pat=dat.dic[[l]]
    
    for(g in 1:length(pt.cuad)){
      
      p.node=p[[g]]
      
      expon=mapply(function(x,y) x^y,p.node,dat.dic.pat,SIMPLIFY = F)
      prod.catgs=lapply(expon,prod)
      denom[g,l]=do.call(prod,prod.catgs)*(w.cuad[g])
      num[g,l]=do.call(prod,prod.catgs)*(w.cuad[g])*pt.cuad[g]
    }
  }
  eap_gpc=colSums(num)/colSums(denom)
  return(eap_gpc)
}

#'@name map_gpc
#'@title Latent Trait Estimation with MAP for gerenalized parcial credit model (gpc) 
#'@description Estimates latent traits for gerenalized parcial credit model (gpc) with MAP
#'data matrix for estimation. Just for unidimensional models.
#'@param estim Output object from the estim_gpc function.
#'@examples
#'\dontrun{
#' data = simulate_polytomous()$data
#' estim=estim_gpc(datos = datos)
#' map_gpc = map_gpc(estim)
#'}
#'@export
map_gpc=function(estim){
	betas=estim$betas
	nitems=length(betas)
	ncatg=estim$ncatg
	dat.dic=estim$dat.dic
	npats=nrow(estim$pats)
	prr=function(thetal,betas,nitems,ncatg){
	  
	  etas=
		lapply(betas,function(x,thetal){x[length(x)]*(thetal-x[1:(length(x)-1)])},thetal=thetal)
	  
	  cumsumetas=lapply(etas,function(x)cumsum(x))
	  den=lapply(cumsumetas,function(x){1+sum(exp(x))})
	  p=list()
	  for(i in 1:nitems){p[[i]]=vector(length = ncatg[i]);
	  p[[i]][1]=1/den[[i]]}
	  
	  for(i in 1:nitems){ 
		for(k in 2:(ncatg[i])){ 
		  p[[i]][k]=exp(etas[[i]][k-1])*p[[i]][k-1]
		}
	  }
	  p<-lapply(p, function(x) ifelse(x <= sqrt(2.2e-16) , sqrt(2.2e-16), x )  ) 
	  p<- lapply(p, function(x) ifelse(x >=  1 - 1e-6,  1 - 1e-6 , x) )
	  
	  # p=lapply(p,function(x){x[-1]})
	  return(p)
	}

	postthetal=function(thetal,l,betas,nitems,ncatg,dat.dic){
	  p=prr(thetal = thetal,betas = betas,nitems = nitems,ncatg = ncatg)
	  pik=mapply(function(x,y)log(x^y),p,dat.dic[[l]],SIMPLIFY = F)  
	  retorno=do.call(sum,lapply(pik,function(x)sum(x)))+log(dnorm(thetal))
	  return(-retorno)
	}
	postthetal(thetal = 0,l = 1,betas = betas,nitems = nitems,ncatg = ncatg,
			   dat.dic=dat.dic)

	map=NULL
	for(l in 1:npats){  
	opt=optim(par = 1,fn = postthetal,l=l,betas=betas,
		   method = "BFGS",
		   nitems=nitems,ncatg=ncatg,dat.dic=dat.dic)
	map[l]=opt$par
	print(l)
	}

	return(map)
}