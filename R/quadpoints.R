#'@name quadpoints
#'@title Quadrature points
#'@export
quadpoints = function(dim, quad_tech = "Gaussian", quad_points = NULL) {
	if ( quad_tech == "QMCEM" ) {

		if ( is.null(quad_points) ) quad_points = 2000

		theta = data.matrix(sobol(n = quad_points, dim = dim, normal = TRUE))
		weights = rep(1, quad_points)
	} else if ( quad_tech == "Gaussian" ) {

		if ( dim == 1 ) quad_points = 40
		else if ( dim == 2 ) quad_points = 20
		else if ( dim == 3 ) quad_points = 10
		else quad_points = 5

		g = gaussHermiteData(quad_points)

		nodes = g$x
		weights = g$w

		lista.nodes = vector("list", length = dim)
		lista.weights = vector("list", length = dim)
		for (i in 1:dim) {
			lista.nodes[[i]] = nodes
			lista.weights[[i]] = weights
		}

		sqrt2 = sqrt(2)
		factor = (sqrt(pi))^dim
		theta = data.matrix( expand.grid(lista.nodes) * sqrt2 )
		weights = apply( expand.grid(lista.weights), MARGIN = 1, function(x) prod(x)) / factor
	}	

	return (list(theta = theta, weights = weights))
}