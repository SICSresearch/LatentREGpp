#include "MulTRI.h"

using namespace Rcpp;

List dichotomous ( IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon,
                   NumericMatrix Rtheta, NumericVector Rweights, IntegerVector Rclusters,
                   NumericMatrix Rinitial_values ) {
  // Converting data types
  irtpp::matrix<char> Y;
  irtpp::matrix<double> theta;
  std::vector<double> weights;
  std::vector<int> clusters;
  irtpp::matrix<double> initial_values;

  irtpp::convert_matrix(Rdata, Y);
  irtpp::convert_matrix(Rtheta, theta);
  irtpp::convert_vector(Rweights, weights);
  irtpp::convert_vector(Rclusters, clusters);
  Rcpp::Rcout << "Ready to convert initial_values\n";
  irtpp::convert_matrix(Rinitial_values, initial_values);

  
  //Estimation object  
  irtpp::dichotomous::estimation e(Y, dim, model, EMepsilon, 
                                   theta, weights, clusters, initial_values);

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

