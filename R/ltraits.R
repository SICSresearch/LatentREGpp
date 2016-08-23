

ltraits = function ( data, dim, model = "2PL", zetas = NULL, 
					quadrature_technique = NULL, quadrature_point = NULL, 
					init_traits = NULL, method = "MAP", by_individuals = TRUE ) {
	# Asserting matrix type
	data = data.matrix(data)

	if ( is.null(zetas) ) {
		print("Item parameters were not loaded\n")
		zetas = lrpp(data, dim, model, quadrature_technique, quadrature_points)$zetas
	} else {
		zeta = data.matrix(zeta)
	}

	q = quadPoints(quadrature_technique, quadrature_points)
	theta = q$theta
	weights = q$weights

	if ( method == "MAP" ) {
		if ( is.null(init_traits) ) 
			traits = ltraitscpp(data = data, dim = dim, model = model, 
								zetas = zetas, theta = theta, weights = weights,
								method = method, by_individuals = by_individuals)
		else {
			traits = ltraitscpp(data = data, dim = dim, model = model, 
								zetas = zetas, theta = theta, weights = weights,
								method = method, by_individuals = by_individuals,
								init_traits = init_traits)
		}
	} else {
		traits = ltraitscpp(data = data, dim = dim, model = model, 
							zetas = zetas, theta = theta, weights = weights,
							method = method, by_individuals = by_individuals)
	}
}