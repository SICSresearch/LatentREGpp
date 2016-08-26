
#'@name standarize
#'@title Standarize function
#'@export
standarize = function( M ) {
	T = matrix(colMeans(M), nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
	C = cov(M - T)
	E = eigen(C)
	V = E$vectors;
	D = diag(E$values)
	S = V%*% sqrt(D)%*%solve(V)
	return (t(solve(S)%*%t(M)))
}