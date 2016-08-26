
#'@name quadpoints
#'@export
quadpoints = function(dim, quadrature_technique = "Gaussian", quadrature_points = NULL) {
	if ( quadrature_technique == "QMCEM" ) {

		if ( is.null(quadrature_points) ) quadrature_points = 2000

		theta = data.matrix(sobol(n = quadrature_points, dim = dim, normal = TRUE))
		weights = rep(1, quadrature_points)
	} else if ( quadrature_technique == "Gaussian" ) {

		if ( dim == 1 ) quadrature_points = 40
		else if ( dim == 2 ) quadrature_points = 20
		else if ( dim == 3 ) quadrature_points = 10
		else quadrature_points = 5

		g = gaussHermiteData(quadrature_points)

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