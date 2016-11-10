#'@name personfit
#'@title Latent traits estimation
#'@description Estimates the latent traits by response pattern according to the Multidimensional Item Response Theory
#'@param data The matrix containing the answers of tested individuals
#'@param dim The dimensionality of the test
#'@param model "1PL", "2PL" or "3PL"
#'@param zetas A matrix with item parameters by item. (can use principal function)
#'@param quad_tech A string with technique. "Gaussian" for Gaussian quadrature 
#'or "QMCEM" for Quasi-Monte Carlo quadrature
#'@param quad_points Amount of quadrature points. If quadratura_technique is "Gaussian". It can be NULL
#'@param init_traits Initial values for MAP latent trait. If EAP it can be NULL
#'@param method "EAP" or "MAP". "MAP" by default.
#'@param by_individuals if True, return latent trait by individual, otherwsie by response pattern. True by default.
#'@param verbose True for get information about estimation process in runtime. False in otherwise.
#'@export
personfit = function ( data, dim, model = "2PL", zetas = NULL, 
					quad_tech = NULL, quad_points = NULL, 
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
	else if ( model == "Bayesian" ) m = 4

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
			return (traits)
		}
	}	
}