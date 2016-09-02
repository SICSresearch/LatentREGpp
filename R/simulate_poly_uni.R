#funcion que simula los datos en poli uni
simulate.poly.uni = function (n = 1000, nitems = 10, ncatgs = c(rep(3, 5), rep(4, 5)), model = "grm", seed_item = 1000,seed_data) 
{
  thetas = list()
  umbrales = list()
  disc=list()
  set.seed(seed_item)
  for ( i in 1:nitems ) {
    beta1 = runif(n = 1, min = -2, max = -1)
    deltas = runif(n = ncatgs[i] - 2, min = -0.4, 1)
    umbrales[[i]] = c(cumsum(c(beta1, exp(deltas))))
    alpha = runif(n = 1, min = 0.7, max = 2)
    thetas[[i]] = c(umbrales[[i]] * alpha, alpha)
    disc[[i]]=alpha
  }
  
  
  
  z <- rnorm(n)#Vector de trazos latentes segun la distribucion
  p <- length(thetas)#numero de items
  nk <- sapply(thetas, length) ##numero de categorias por item
  
  prob <- if (model == "grm") { #aca si entra
    thetas <- lapply(thetas, function(x) {
      n_x <- length(x)
      cbind(plogis(matrix(x[-n_x], n, n_x - 1, TRUE) - x[n_x] * z), 1)
    })
    lapply(thetas, function(x) {
      nc <- ncol(x)
      cbind(x[, 1], x[, 2:nc] - x[, 1:(nc - 1)])
    })
  }
  
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
  
  retorno = list(data = data, params.it = mapply(function(x,y){c(x,y)},disc,d))
  return(retorno)
}
