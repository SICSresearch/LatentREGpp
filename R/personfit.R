#'@name personfit
#'@title Latent traits estimation
#'@description Estimates the latent traits by using either the Expected A Posteriori (EAP) 
#'or Mode A Posteriori (MAP) method. A Normal distribution with mean vector zero and 
#'covariance matrix the identity is assumed. Quasi-Monte Carlo quadrature is advised
#'when the data dimension is large \eqn{(>3)}.
#'@param data The matrix containing the answers of tested individuals
#'@param dim The dimensionality of the test
#'@param model "1PL", "2PL" or "3PL"
#'@param zetas The item parameters. A matrix of dim (num of items * num of parameters
#' from item that has the greater number of categories) where each row is a vector
#' of the form: 
#' 
#' \eqn{(\boldsymbol{a}_i,\gamma_{i1},\gamma_{i2},...,\gamma_{im_i},
#' NA,NA,...,NA,c_i)} according to the notation in the section "Notation". 
#' The function "LatentREGpp::itemfit( )" returns the zetas with this structure.
#'@param quad_tech A string with technique. "Gaussian" for Gaussian quadrature 
#'or "QMCEM" for Quasi-Monte Carlo quadrature
#'@param quad_points Amount of quadrature points by dimention. Default NULL
#'@param by_individuals if True, return latent trait by individual, otherwsie by response pattern. True by default.
#'@param init_traits Initial values by pattern or by individual as the case may be. Default NULL. 
#'@param method "EAP" or "MAP". "MAP" by default.
#'@param verbose True for get information about estimation process in runtime. False in otherwise.
#'
#'@section Methods to estimate the Latent Traits:
#' \describe{
#'   \item{EAP}{
#'   In general the EAP Method is based on the next expression. 
#'   \deqn{\frac{\int \theta_lp(U_{l}=u_l /\theta_l,\zeta)
#'   p(\theta_l/\eta)\partial\theta_l}{\int p(U_{l}=u_l /\theta_l,\zeta)
#'   p(\theta_l/\eta)\partial\theta_l}}
#'   where:
#'   \describe{
#'   \item{} \eqn{\theta_l} is the latent trait associated with pattern l.
#'   \item{} \eqn{U_l} refers to response pattern l
#'   \item{} \eqn{\zeta} are items parameters
#'   \item{} \eqn{\eta} are the hiperparameters from the prior distribution to the traits. }  }
#'   
#'   \item{MAP} {The method consists of maximize the following expression regard to
#'   \eqn{\theta_l}
#'   \deqn{\frac{p(U_{l}=u_l /\theta_l,\zeta)p(\theta_l/\eta)} {\int p(U_{l}=u_l /\theta_l,\zeta)
#'   p(\theta_l/\eta)\partial\theta_l}}}
#'
#'}
#'
#'@section Notation:
#'
#'In the Polytomous Multidimensional case, the probability that an examinee with latent trait vector 
#'\eqn{\theta_l} responses categorie k to item i is,
#'\deqn{ P(U_{li}=k \mid \boldsymbol{\theta}_l,
#'\boldsymbol{a}_i,\gamma_{ik}) =c_i+(1-c_i)\boldsymbol{\Psi}(\eta_{lik})}
#'where 
#'\describe{
#'\item{} {\eqn{\eta_{lik}=\boldsymbol{a}_i^t\boldsymbol{\theta}_l + \gamma_{ik}}}
#'}
#'\describe{
#'\item{}{\eqn{\boldsymbol{a}_i=(\boldsymbol{a}_1,\boldsymbol{a}_2,...,\boldsymbol{a}_d)^t} 
#'is a parameter associated with discrimination of item, d is the data dimention }
#'\item{}{\eqn{\boldsymbol{\theta}_l=(\boldsymbol{\theta}_{l1},\boldsymbol{\theta}_{l2},...,\boldsymbol{\theta}_{ld})^t}}
#'is the latent trait multidimensional associated with pattern l.
#'\item{}{\eqn{\gamma_{ik}} is a parameter associated with item and the categorie,
#' \eqn{i=1,2,..,p} (number of items) \eqn{k=1,2...,m_i} (number of categories from item)}
#' \item{}{\eqn{c_i} is the guessing parameter, make sense in the dichotomous case.}
#'}
#'
#'@return depends on value of by_individuals argument.
#'\describe{
#'\item{}{If by_individuals=TRUE. Returns a matrix with latent trait for each individual}
#'\item{}{If by_individuals=FALE. Returns a list of length 3 with the latent traits for each pattern,
#'the reponse patterns and the frecuencie of each pattern}
#'}
#' @references
#' lllll
#' 
#' llll
#' @examples
#' 
#' #Example 1
#' 
#' #simulate 10 polyotmous items, the first 5 with 4 response categories and the
#' #others with 5 response categories, By default the number of individuals is 1000
#' #, model is "2PL" and the dimention is 1
#' dats=simulate_polytomous(ncatgs = c(rep(4,10),rep(5,10)),seed_data = 5000L)
#' #estimate the items parameters
#' est=itemfit(dats$data,dim = 1,model = "2PL")
#' est$zetas
#' #calculates the latent traits estimation.
#' personfit(dats$data,dim = 1,zetas = est$zetas,method = "MAP")
#' 
#' #Example 2
#' 
#' #simulate 10 dichotomous items, the trait dimention is 2
#' #other arguments by default
#' dats=simulate_dichotomous(dim.data = 2,seed_data = 5000L)
#' #estimate the items parameters
#' est=itemfit(dats$data,dim = 2,model = "2PL")
#' est$zetas
#' #calculates the latent traits estimation
#' personfit(dats$data,dim = 2,zetas = est$zetas,method = "MAP",by_individuals=F)
#'@export
personfit = function ( data, dim, model = "2PL", zetas , 
                       quad_tech = "Gaussian", quad_points =NULL, 
                       init_traits = NULL, method = "MAP", by_individuals = TRUE,
                       verbose = FALSE ) {
  # Quadrature technique
  if ( is.null(quad_tech) ) {
    if ( dim < 5 ) quad_tech = "Gaussian"
    else quad_tech = "QMCEM"
  } else {
    if ( quad_tech != "Gaussian" && quad_tech != "QMCEM" )
      stop("Quadrature technique not found")
    
    if ( dim >= 5 && quad_tech == "Gaussian" ) {
      message("For dim >= 5 QMCEM quadrature technique is recommended")
      input = readline(prompt = "Are you sure you want continue with Gaussian quadrature? [y/n]: ")
      if ( input == "n" || input == "N" )
        quad_tech = "QMCEM"
      else quad_tech = quad_tech
    }
    else quad_tech = quad_tech
  }
  
  # Asserting matrix type
  data = data.matrix(data)
  
  # Type of data
  dichotomous_data = is_data_dicho(data)
  if ( dichotomous_data == -1 )
    stop("Inconsistent data")
  
  # Model id
  if ( model == "1PL" ) m = 1
  else if ( model == "2PL" ) m = 2
  else if ( model == "3PL" ) m = 3
  
  q = quadpoints(dim = dim, quad_tech = quad_tech, 
                 quad_points = quad_points)
  theta = q$theta
  weights = q$weights
  
  if ( is.null(zetas) ) 
    zetas = itemfit(data = data, dim = dim, model = model,
                    quad_tech = quad_tech, quad_points = quad_points,
                    verbose = verbose)$zetas
  
  if ( method == "MAP" && !is.null(init_traits) ) {
    init_traits = data.matrix(init_traits)
    traits = personfitcpp(Rdata = data, dim = dim, model = m, 
                          Rzetas = zetas, Rtheta = theta, Rweights = weights,
                          method = method, by_individuals = by_individuals,
                          dichotomous_data = dichotomous_data,
                          Rinit_traits = init_traits)
    
  } else {
    if ( method == "MAP" ) {
      traits = personfitcpp(Rdata = data, dim = dim, model = m, 
                            Rzetas = zetas, Rtheta = theta, Rweights = weights,
                            method = "EAP", by_individuals = FALSE,
                            dichotomous_data = dichotomous_data)
      traits$latent_traits = standarize(traits$latent_traits)
      traits = personfitcpp(Rdata = data, dim = dim, model = m, 
                            Rzetas = zetas, Rtheta = theta, Rweights = weights,
                            method = method, by_individuals = by_individuals,
                            dichotomous_data = dichotomous_data,
                            Rinit_traits = traits$latent_traits)
    } else {
      traits = personfitcpp(Rdata = data, dim = dim, model = m, 
                            Rzetas = zetas, Rtheta = theta, Rweights = weights,
                            method = method, by_individuals = by_individuals,
                            dichotomous_data = dichotomous_data)
      traits$latent_traits = standarize(traits$latent_traits)
    }
  }
  return (traits)
}