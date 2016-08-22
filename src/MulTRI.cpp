#include <Rcpp.h>
#include "polytomous/estimation/estimation.h"
#include "dichotomous/estimation/estimation.h"
#include "util/constants.h"
#include "util/matrix.h"
#include "util/general.h"

using namespace Rcpp;

//' Estimates the parameters of a test
//'
//' @param RData Input dataset.
//' @export
// [[Rcpp::export]]
List dichotomous ( IntegerMatrix RData, unsigned int dim, int model, double EMepsilon,
                              std::string wd ) {
  irtpp::matrix<char> Y;
  convert_matrix(RData, Y);
  irtpp::setwd(wd);
  irtpp::dichotomous::estimation e(Y, dim, model, EMepsilon);

  //convert_matrix(quad, e.estimation_data.theta);
  // convert_matrix(weights, e.estimation_data.weights);
  // e.estimation_data.G = G;
  // e.build_matrixes();
  e.EMAlgorithm();

  NumericMatrix zetas(e.data.p, e.data.d + 2);

  int current_zeta = e.get_iterations() % 3;

  for ( int i = 0; i < e.data.p; ++i ) {
    int j = 0;
    if ( e.data.m.parameters > 1 )
      for ( ; j < e.data.d; ++j ) zetas(i,j) = e.data.zeta[current_zeta][i](j);
    zetas(i,j) = e.data.zeta[current_zeta][i](j);
    if ( e.data.m.parameters == 3 ) zetas(i,j+1) = e.data.zeta[current_zeta][i](j+1);
    else zetas(i,j+1) = 0;
  }
  double loglik = e.log_likelihood();
  List output = List::create(zetas,loglik);
  return output;
}


//' Estimates the parameters of a test
//'
//' @param RData Input dataset.
//' @export
// [[Rcpp::export]]
List polytomous ( IntegerMatrix RData, unsigned int dim, int model, double EMepsilon,
                              std::string wd ) {
  irtpp::matrix<char> Y;
  convert_matrix(RData, Y);
  irtpp::setwd(wd);
  irtpp::dichotomous::estimation e(Y,dim);

  //convert_matrix(quad, e.estimation_data.theta);
  // convert_matrix(weights, e.estimation_data.weights);
  // e.estimation_data.G = G;
  // e.build_matrixes();
  e.EMAlgorithm();

  NumericMatrix zetas(e.data.p,e.data.d+2);

  int current_zeta = e.get_iterations() % 3;

  for ( int i = 0; i < e.data.p; ++i ) {
    int j = 0;
    if ( e.data.m.parameters > 1 )
      for ( ; j < e.data.d; ++j ) zetas(i,j) = e.data.zeta[current_zeta][i](j);
      zetas(i,j) = e.data.zeta[current_zeta][i](j);
      if ( e.data.m.parameters == 3 ) zetas(i,j+1) = e.data.zeta[current_zeta][i](j+1);
      else zetas(i,j+1) = 0;
  }
  double loglik = e.log_likelihood();
  List output = List::create(zetas,loglik);
  return output;
}
