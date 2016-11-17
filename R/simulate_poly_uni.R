#funcion que simula los datos en poli uni
simulate.poly.uni = function (n=1000, nitems=20, ncatgs=rep(2,20), seed_item=5000L,seed_data=2000L,model="2PL") 
{
  thetas = list()
  umbrales = list()
  disc=list()
  set.seed(seed_item)
  for ( i in 1:nitems ) {
    beta1 = runif(n = 1, min = -1.4, max = -0.2)
    deltas = runif(n = ncatgs[i] - 2, min = 0,max = 1)
    umbrales[[i]] = c(cumsum(c(beta1, exp(deltas))))
    if(model!="1PL"){alpha = runif(n = 1, min = 0.7, max = 2)}
    else{alpha=1}
    thetas[[i]] = c(umbrales[[i]] * alpha, alpha)
    disc[[i]]=alpha
  }
  
  z <- rnorm(n)#Vector de trazos latentes segun la distribucion
  p <- length(thetas)#numero de items
  nk <- sapply(thetas, length) ##numero de categorias por item
  
  thetas <- lapply(thetas, function(x) {
    n_x <- length(x)
    cbind(plogis(matrix(x[-n_x], n, n_x - 1, TRUE) - x[n_x] * z), 1)
  })
  
  prob <-
    lapply(thetas, function(x) {
      nc <- ncol(x)
      cbind(x[, 1], x[, 2:nc] - x[, 1:(nc - 1)])
    })
  
  if(model=="3PL"){
    set.seed(2000L)
    cc=runif(p,0,0.25)
    for(i in 1:p){
      prob[[i]][,2]=cc[i]+((1-cc[i])*prob[[i]][,2])
      prob[[i]][,1]=1-prob[[i]][,2]
    }
  }
  else{cc=rep(0,p)}
  
  data <- matrix(0, n, p)
  set.seed(seed_data)
  for (j in 1:p) {  
    ##Extrae las muestras con las probabilidades halladas anteriormente segun el modelo
    for (i in 1:n) data[i, j] <- sample(nk[j], 1, prob = prob[[j]][i, ])
  }
  
  d = umbrales
  for ( j in 1:p ) {
    d[[j]] = -1 * d[[j]] * disc[[j]]
  }
  
  itempars=mapply(function(x,y){c("disc"=x,"-dif*disc"=y)},disc,d)
  retorno = list(data = data, params.it = itempars,traits=z,guessing=cc)
  return(retorno)
}
