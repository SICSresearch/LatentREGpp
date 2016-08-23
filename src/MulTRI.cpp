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
                   NumericMatrix theta, NumericVector weights, IntegerVector Rclusters,
                   std::string initial_values ) {
  irtpp::matrix<char> Y;
  irtpp::convert_matrix(Rdata, Y);
  std::vector<int> clusters;
  irtpp::convert_vector(Rclusters, clusters);

  //Estimation config
  if(initial_values.empty()) {
    initial_values = irtpp::NONE;
    //std::cout<<"Enter here"<<std::endl;
    //std::cout<<initial_values<<std::endl;
  }
  
  irtpp::dichotomous::estimation e(Y, dim, model, EMepsilon, clusters,
                                   irtpp::GAUSSIAN, irtpp::DEFAULT_QMCEM_POINTS, 
                                   irtpp::EMPTY_INTEGER_VECTOR, initial_values);
  
  //Quadrature points config 
  irtpp::convert_matrix(theta, e.data.theta);
  irtpp::convert_vector(weights, e.data.w);
  e.data.G = e.data.theta.rows();
  e.build_matrixes();

  //EM
  e.EMAlgorithm();

  NumericMatrix zetas(e.data.p, e.data.d + 2);
  int current_zeta = e.get_iterations() % irtpp::ACCELERATION_PERIOD;

  for ( int i = 0; i < e.data.p; ++i ) {
    int j = 0;
    if ( e.data.m.parameters > 1 )
      for ( ; j < e.data.d; ++j ) zetas(i,j) = e.data.zeta[current_zeta][i](j);
    zetas(i,j) = e.data.zeta[current_zeta][i](j);
    if ( e.data.m.parameters == 3 ) zetas(i,j+1) = e.data.zeta[current_zeta][i](j+1);
    else zetas(i,j+1) = 0;
  }

  double loglikelihood = e.log_likelihood();
  List output = List::create(Rcpp::Named("zetas") = zetas,
                             Rcpp::Named("Loglikelihood") = loglikelihood);

  return output;
}

