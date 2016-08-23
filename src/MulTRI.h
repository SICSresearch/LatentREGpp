#include <Rcpp.h>
#include "polytomous/estimation/estimation.h"
#include "dichotomous/estimation/estimation.h"
#include "util/constants.h"
#include "util/matrix.h"
#include "util/general.h"

using namespace Rcpp;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::plugins(openmp)]]

// Estimates the parameters of a test
//
// @param Rdata Input dataset.
// @param dim Model Dimension
// @param model 1 2 or 3PL
// @param EMepsilon Convergence value for the algorithm
// @param theta Quadrature Points
// @param weights Quadrature Points Weights
// [[Rcpp::export]]
List dichotomous ( IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon,
                   NumericMatrix Rtheta, NumericVector Rweights, 
                   IntegerVector Rindividual_weights,
                   IntegerVector Rclusters = IntegerVector::create(),
                   NumericMatrix Rinitial_values = NumericMatrix(0, 0, 0) );