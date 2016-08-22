#######################################################################
#' @name MulTRI
#' @docType package
#' @title MulTRI : Item Response Theory Implemented in R and Cpp
#' @description Mul-Tri is a c++ implementation of the Multidimensional Item Respone Theory (MIRT)
#' cappable of performing parameter and traits estimations. It also provides a list of options to 
#' perform optiman analysis and provides usefull information about the obtained model.
#' @details
#' \tabular{ll}{
#'Package: \tab MulTRI\cr
#'Type: \tab Package\cr
#'Version: \tab 0.0.5\cr
#'Date: \tab 2016-08-22\cr
#'License: \tab MIT + file LICENSE \cr
#'}
#'@author SICS Research Team
#'@keywords IRT MIRT Psychometry 
#'@useDynLib MulTRI
#'@importFrom Rcpp sourceCpp
#'@importFrom randtoolbox sobol
#'@importFrom fastGHQuad gaussHermiteData
#'@section Getting Started:
#'Get started with the MulTRI package browsing the index of this documentation
#'if you need help the vignettes should be helpful.
#'@section Getting Started:
#'The IRTpp package allows you to use the MulTRI methodology for simulating, analyzing and scoring tests \cr
#'You can browse the package vignettes to get started.
#'
NULL

#'@name multri
#'@title Parameter estimation of a test
#'@description Estimates the test parameters according to the Multidimensional Item Response Theory
#'@param data The matrix containing the answers of tested individuals
#'@param dim The Dimensionality of the test
#'@param model 1PL, 2PL or 3PL
#'@param EMepsilon Convergence value to determine the accuracy of the test
#'@param clusters Clusters per dimention
#'@param quadratura_technique Quasi-Monte Carlo or Gaussian
#'@param quadrature_points Amount of quadrature points
#'@param individual_weights Weights of the quadrature points
#'@param initial_values Initial Values of the estimation
#'@export
multri = function(data, dim, model = "2PL", EMepsilon = 0.0001, clusters = NULL,
				  quadrature_technique = NULL, quadrature_points = NULL, 
				  indivudual_weights = NULL,
				  initial_values = NULL ) {
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

		theta = data.matrix(sobol(n = quadrature_points, dim = dim, normal = TRUE))
		weights = rep(1, quadrature_points)
	} else if ( quadrature_technique == "Gaussian" ) {

		if ( dim == 1 ) quadrature_points = 40
		else if ( dim == 2 ) quadrature_points = 20
		else if ( dim == 3 ) quadrature_points = 10
		else quadrature_points = 5

		g = gaussHermiteData(quadrature_points)

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
		# TODO Sergio's initial values here
		# TODO Find pinned items

		# Item parameters estimation
		if ( dichotomous_data )
			dichotomous(RData = data, dim = dim, model = m, EMepsilon = EMepsilon, 
						theta = theta, weights = weights)
		#else
		#	poly
	}
}