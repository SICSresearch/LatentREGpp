# Main file to call R functions
# Package main function

MulTRI = function(data, dim, model = "2PL", EMepsilon = 0.0001, clusters = NULL,
				  quadrature_technique = NULL, quadrature_points = NULL, 
				  indivudual_weights = NULL,
				  starting_values = NULL ) {
	# Asserting matrix type of data
	data = data.matrix(data)

	# Model id
	m = 2
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

	# Quadrature points and weights
	theta = NULL
	weights = NULL
	if ( quadrature_technique == "QMCEM" ) {

		if ( is.null(quadrature_points) ) quadrature_points = 2000

		theta = data.matrix(randtoolbox::sobol(n = quadrature_points, dim = dim, normal = TRUE))
		weights = rep(1, quadrature_points)
	} else if ( quadrature_technique == "Gaussian" ) {

		if ( dim == 1 ) quadrature_points = 40
		else if ( dim == 2 ) quadrature_points = 20
		else if ( dim == 3 ) quadrature_points = 10
		else quadrature_points = 5

		g = fastGHQuad::gaussHermiteData(quadrature_points)

		nodes = g$x
		weigths = g$w

		lista.nodes = vector("list", length = dim)
		lista.weigths = vector("list", length = dim)
		for (i in 1:dim) {
			lista.nodes[[i]] = nodes
			lista.weigths[[i]] = weigths
		}

		sqrt2 = sqrt(2)
		factor = (sqrt(pi))^dim
		theta = data.matrix( expand.grid(lista.nodes) * sqrt2 )
		weights = apply( expand.grid(lista.weigths), MARGIN = 1, function(x) prod(x)) / factor
	}


	# Type of data
	dichotomous_data = TRUE
	# Find here the data type
	# dichotomous_data = data_type(data)


	if ( dim == 1 ) {
		# Item parameters estimation
		if ( dichotomous_data )
			dichotomous(RData = data, dim = dim, model = m, EMepsilon = EMepsilon, 
						theta = theta, weights = weights)
		#else
		#	poly
	} else {
		# TODO Sergio's starting values here
		# TODO Find pinned items

		# Item parameters estimation
		if ( dichotomous_data )
			dichotomous(RData = data, dim = dim, model = m, EMepsilon = EMepsilon, 
						theta = theta, weights = weights)
		#else
		#	poly
	}
}