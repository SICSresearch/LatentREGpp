#'@name ltraits
#'@title Latent traits estimation
#'@export
ltraits = function ( data, dim, model = "2PL", zetas = NULL, 
					quadrature_technique = NULL, quadrature_points = NULL, 
					init_traits = NULL, method = "MAP", by_individuals = TRUE ) {
	# Asserting matrix type
	data = data.matrix(data)

	if ( is.null(zetas) ) {
		print("Item parameters were not loaded\n")
		zetas = latentreg(data = data, dim = dim, model = model, quadrature_technique = quadrature_technique, 
			quadrature_points = quadrature_points)$zetas
	} else {
		zetas = data.matrix(zetas)
	}

	# Model id
	if ( model == "1PL" ) m = 1
	else if ( model == "2PL" ) m = 2
	else if ( model == "3PL" ) m = 3

	# Quadrature technique
	if ( is.null(quadrature_technique) ) {
		if ( dim < 5 ) quadrature_technique = "Gaussian"
		else quadrature_technique = "QMCEM"
	} else {
		if ( dim >= 5 && quadrature_technique == "Gaussian" ) {
			print("Better use QMCEM")
			# Try to change the qudrature technique
		}
	}

	q = quadpoints(dim = dim, quadrature_technique = quadrature_technique, 
				   quadrature_points = quadrature_points)
	theta = q$theta
	weights = q$weights

	# Type of data
	dichotomous_data = is_data_dicho(data)
	if ( dichotomous_data == -1 ) {
		print("Inconsistent data")
		return -1
	}

	if ( method == "MAP" && !is.null(init_traits) ) {
		init_traits = data.matrix(init_traits)
		traits = ltraitscpp(Rdata = data, dim = dim, model = m, 
								Rzetas = zetas, Rtheta = theta, Rweights = weights,
								method = method, by_individuals = by_individuals,
								dichotomous_data = dichotomous_data,
								Rinit_traits = init_traits)

	} else {
		if ( method == "MAP" ) {
			traits = ltraitscpp(Rdata = data, dim = dim, model = m, 
								Rzetas = zetas, Rtheta = theta, Rweights = weights,
								method = "EAP", by_individuals = FALSE,
								dichotomous_data = dichotomous_data)
			traits$latent_traits = standarize(traits$latent_traits)
			traits = ltraitscpp(Rdata = data, dim = dim, model = m, 
								Rzetas = zetas, Rtheta = theta, Rweights = weights,
								method = method, by_individuals = by_individuals,
								dichotomous_data = dichotomous_data,
								Rinit_traits = traits$latent_traits)
		} else {
			traits = ltraitscpp(Rdata = data, dim = dim, model = m, 
								Rzetas = zetas, Rtheta = theta, Rweights = weights,
								method = method, by_individuals = by_individuals,
								dichotomous_data = dichotomous_data)
			traits$latent_traits = standarize(traits$latent_traits)
			return (traits)
		}
	}	
}