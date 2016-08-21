# Main file to call R functions
# Package main function

MulTRI = function(data, dim, model = "2PL", EMepsilon = 0.0001, clusters = NULL,
				  quadrature_technique = "Gaussian", weights = NULL,
				  starting_values = NULL ) {
	wd = paste( c(system.file(package = "IRTPP"), "/"), collapse = "" )
	MulTRICall(data, dim, wd)
}