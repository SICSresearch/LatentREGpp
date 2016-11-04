#'@name quadpoints
#'@title Quadrature points
#'@description Return a list with quadrature points according dimensionality, technique
#'and number of points.
#'@param dim Dimension of the quadrature
#'@param quad_tech A string with technique. "Gaussian" for Gaussian quadrature. 
#'or "QMCEM" for Quasi-Monte Carlo quadrature.
#'@param quad_points Amount of quadrature points. If quadratura_technique is "Gaussian". 
#'It can be NULL, in Quasi-Monte Carlo it is 2000 by default. 
#'@examples
#'\dontrun{qp = quadpoints(dim = 4,quad_tech = "QMCEM",quad_points = 3000)}
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