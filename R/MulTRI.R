# Main file to call R functions
# Package main function

MulTRI = function(data, dim, model = "2PL", EMepsilon = 0.0001, clusters = NULL,
				  quadrature_technique = "Gaussian", weights = NULL,
				  starting_values = NULL ) {
	# Working directory
	wd = paste( c(system.file(package = "IRTPP"), "/"), collapse = "" )

	# Model id
	m = 2
	if ( model == "1PL" ) m = 1
	else if ( model == "2PL" ) m = 2
	else if ( model == "3PL" ) m = 3

	# Type of data
	dichotomous_data = TRUE
	#dichotomous_data = data_type(data)

	if ( dim == 1 ) {
		# Item parameters estimation
		if ( dichotomous_data )
			dichotomous(RData = data, dim = dim, model = m, EMepsilon = EMepsilon, 
								wd = wd)
		else 
			polytomous(RData = data, dim = dim, model = m, EMepsilon = EMepsilon, 
								wd = wd)
	} else {
		# FInd starting values here

		# Item parameters estimation
		if ( dichotomous_data )
			dichotomous(RData = data, dim = dim, model = m, EMepsilon = EMepsilon, 
								wd = wd)
		else 
			polytomous(RData = data, dim = dim, model = m, EMepsilon = EMepsilon, 
								wd = wd)

	}
}