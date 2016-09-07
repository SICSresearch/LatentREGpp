#'@name standarize
#'@title Standarize function
#'@description This function takes as input any data matrix and standarizes it.
#'@param M A matrix with data
#'@export
standarize = function( M ) {
	d = ncol(M)
	T = matrix(colMeans(M), nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
	C = cov(M)
	E = eigen(C)
	V = E$vectors
	if ( d > 1 )
		D = diag(E$values)
	else
		D = E$values
	S = V%*%sqrt(D)%*%t(V)
	return (t(solve(S)%*%t(M - T)))
}