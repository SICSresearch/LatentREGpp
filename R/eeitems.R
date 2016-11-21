probqi=function(betasi,nodes,item,ncatg,dicomod=dicomod){
  
  ind1 <- c(1, cumsum(ncatg[-length(ncatg)]) + 1)
  ind2 <- cumsum(ncatg)
  dimdata=ncol(nodes)

  
  pr=matrix(NA,nrow=nrow(nodes),ncol=length(c(ind1[item]:ind2[item])))
  for(i in 1:nrow(nodes)){  
    alphai=betasi[1:dimdata]
    gamma=betasi[(dimdata+1):length(betasi)]
    if(dicomod==T){gamma=betasi[(dimdata+1):(length(betasi)-1)]}
    eta=alphai%*%nodes[i,]+gamma
    x=1/(1+exp(-eta))
    past=c(1,x,0)  
    pr[i,]=-diff(past)
    if(dicomod==T){
      cc=betasi[length(betasi)]
      pr[i,]=cc+((1-cc)*pr[i,])
      }
  }
  
  pr<-ifelse(pr <= sqrt(2.2e-16) , sqrt(2.2e-16), pr) 
  pr<- ifelse(pr >=  1 - 1e-16,  1 - 1e-16,pr )
  
  return(pr)
  
}

ssdd=function(betas,r,nodes,ncatg,f=NULL,dicomod=FALSE){
  nitems=nrow(betas)
  
  Qi=function(betasi,ritem,nodes,item,ncatg,dicomod){
    prr=probqi(betasi = betasi,nodes = nodes,item = item,ncatg = ncatg,dicomod=dicomod)
    sum(prr*ritem)
  }

  Qidico=function(betasi,nodes,item,ncatg,ritem,fitem,dicomod){
    prr=probqi(betasi = betasi,nodes = nodes,item = item,
      ncatg = ncatg,dicomod=dicomod)[,2]
    sum(ritem*log(prr)+(fitem-ritem)*log(1-prr))
  }

  inf=list()
  for(item in 1:nitems){  #item=1
    if(dicomod){
      betasi=na.omit(betas[item,])
      ritem=r[,item]
      # opt=optim(par = betasi,fn = Qidico,ritem=ritem,item=item,
      #            ncatg=ncatg,nodes=nodes,fitem=f,dicomod=dicomod,method = "CG")
      Hess=hessian(func = Qidico,x = betasi,ritem=ritem,item=item,
              ncatg=ncatg,nodes=nodes,fitem=f,dicomod=dicomod)
      
      if(sum((eigen(Hess)$values>=0))!=0){
        print("The information matrix is not positive definite for the item:")
        print(item)
      }
      
      inf[[item]]=solve(-Hess)
    }
    else{
      betasi=na.omit(betas[item,])
      ritem=do.call(rbind,lapply(r,function(x){na.omit(x[item,])}))
      # opt=optim(par = betasi,fn = Qi,ritem=ritem,item=item,
      #           ncatg=ncatg,nodes=nodes,dicomod=dicomod)
      Hess=hessian(func = Qi,x = betasi,ritem=ritem,item=item,
              ncatg=ncatg,nodes=nodes,dicomod=dicomod)
      
      if(sum((eigen(Hess)$values>=0))!=0){
        # print(eigen(Hess)$values)
        print("The information matrix is not positive definite for the item:")
        print(item)
        }
      
      inf[[item]]=solve(-Hess)
    } 
   }
  return(inf)
}

