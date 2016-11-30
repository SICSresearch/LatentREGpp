#'@name itemfit_bayesian
#'@title Bayesian parameter estimation of a test
#'@description Estimates the test parameters according to the Multidimensional Item Response Theory with
#'bayesian adjust for dichotomous data
#'@param data The matrix containing the answers of tested individuals
#'@param dim The dimensionality of the test
#'@param model "1PL", "2PL" or "3PL"
#'@param EMepsilon Convergence value to determine the accuracy of the test
#'@param clusters A vector with cluster per dimension 
#'@param quad_tech A string with technique. "Gaussian" for Gaussian quadrature 
#'or "QMCEM" for Quasi-Monte Carlo quadrature. If NULL it's selected according to the model's dimension (QMCEM if dim>3).
#'@param quad_points Amount of quadrature points. If quadratura_technique is "Gaussian". It can be NULL
#'@param individual_weights A vector with Weights of the quadrature points.
#'@param initial_values A matrix with initial values for estimation process. Be sure about
#'dimension, model and consistency with data. 
#'@param noguessing In 3PL model and dimension is greater than 1, If true, guessing parameter will not be estimated in zeta vector. Instead
#'c value will have a default initial value. Otherwise guessing parameter will be estimated with zeta vector.
#'@param verbose True for get information about estimation process in runtime. False in otherwise. 
#'@param save_time True for save estimation time. False otherwise.
#'@section Model:
#'\describe{
#'  Bayesian model is based in itemfit models. It has a \eqn{Q_{i}} function to optimize according parameters like in itemfit. However this model 
#'  is given by: 
#'  \deqn{Q_{i} = N * log(P_{\zeta_{i}}(\zeta_{i})) + \hat{Q_{i}}} 
#'
#'  Where i index is referenced for items in test.
#'
#'  Then, log posterior is given by:
#'
#'  \deqn{log(P_{\zeta_{i}}(\zeta_{i})) = - \frac{N}{2} (\frac{(a_{1i} - \mu_{a1i})^2}{\sigma_{1i}^2} + \cdots + \frac{(a_{Di} - \mu_{aDi})^2}{\sigma_{Di}^2} + \frac{(d_{i} - \mu_{di})^2}{\sigma_{di}^2} + \frac{(c_{i} - \mu_{ci})^2}{\sigma_{ci}^2}) }
#'
#'  Where a,d and c are parameters, D is the dimension of test. You can give the \eqn{\mu} values for each parameters through initial values matrix. In otherwise \eqn{\mu} will have default initial values value
#'  \eqn{\sigma^2} values are constant \eqn{\sigma_a^2 = 0.64}, \eqn{\sigma_d^2 = 4}, \eqn{\sigma_c^2 = 0.009}
#'}
#'@export
itemfit_bayesian = function(data, dim, model = "2PL", EMepsilon = 1e-4, clusters = NULL,
				  quad_tech = NULL, quad_points = NULL, 
				  individual_weights = as.numeric(c()),
				  initial_values = NULL, noguessing = TRUE, 
				  verbose = TRUE, save_time = TRUE ) {
	
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
		}
	}

	if ( save_time ) first_time = Sys.time()

	# Asserting matrix type of data
	data = data.matrix(data)

	# Type of data
	dichotomous_data = is_data_dicho(data)
	if ( dichotomous_data == -1 )
		stop("Inconsistent data")

	# Model id
	#in this function model always has a 4 value
	#but type has a model value
	#TODO
	#m = 4

	if ( model == "1PL" ) m = 1
	else if ( model == "2PL" ) m = 2
	else if ( model == "3PL" ) m = 3
	else stop("Model not valid")

	q = quadpoints(dim = dim, quad_tech = quad_tech, 
				   quad_points = quad_points)
	theta = q$theta
	weights = q$weights

	if ( dim == 1 ) {
		# Item parameters estimation
		obj_return = itemfitcpp_bayesian(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
							Rtheta = theta, Rweights = weights, 
							Rindividual_weights = individual_weights,
							dichotomous_data = dichotomous_data, 
							verbose = verbose )
	} else {
		if ( dichotomous_data ) {

			if ( is.null(clusters) ) {
				#1. Find a temporal cluster
				p = ncol(data)
				d = dim
				clusters = find_temporal_cluster(p=p,d=d)
			} else
				if ( length(clusters) != dim )
					stop("Clusters length must be equal to the number of dimensions")

			#Initial values
			if ( is.null(initial_values) ) {
				# Item parameters estimation with no initial values provided
				obj_return = itemfitcpp_bayesian(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
									Rtheta = theta, Rweights = weights, 
									Rindividual_weights = individual_weights,
									dichotomous_data = dichotomous_data,
									Rclusters = clusters,
									noguessing = noguessing,
									verbose = verbose)
			} else {
				initial_values = data.matrix(initial_values)
				if ( nrow(initial_values) != ncol(data) )
					stop("Inconsistent initial_values. Number of rows must be equal to the number of items")

				# Item parameters estimation 
				obj_return = itemfitcpp_bayesian(Rdata = data, dim = dim, model = m, EMepsilon = EMepsilon, 
									Rtheta = theta, Rweights = weights, 
									Rindividual_weights = individual_weights,
									dichotomous_data = dichotomous_data,
									Rclusters = clusters,
									Rinitial_values = initial_values, 
									noguessing = noguessing,
									verbose = verbose)
			}
		  
		} else {
			stop("Bayesian model is not available for polytomous data")
		}
	}

	if ( save_time ) {
		second_time = Sys.time()
		obj_return$time = second_time - first_time	
	}

	obj_return$dimension = dim
	obj_return$model = model
	obj_return$clusters = clusters
	obj_return$convergence = obj_return$iterations < 500
	obj_return$epsilon = EMepsilon
	obj_return$quadpoints = q
  	
  	if ( dichotomous_data )
		colnames(obj_return$zetas) = c(paste("a",c(1:obj_return$dimension),sep = ""),"d","c")
	else {
		stop("You can not have an estimation in bayesian with polytomous data")
	}
  
	return (obj_return)
}