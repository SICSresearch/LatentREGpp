#include "MulTRI.h"

using namespace Rcpp;

List dichotomous ( IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon,
                   NumericMatrix Rtheta, NumericVector Rweights, 
                   IntegerVector Rindividual_weights, IntegerVector Rclusters,
                   NumericMatrix Rinitial_values ) {
  // Converting data types
  irtpp::matrix<char> Y;
  irtpp::matrix<double> theta;
  std::vector<double> weights;
  std::vector<int> individual_weights;
  std::vector<int> clusters;
  irtpp::matrix<double> initial_values;

  irtpp::convert_matrix(Rdata, Y);
  irtpp::convert_matrix(Rtheta, theta);
  irtpp::convert_vector(Rweights, weights);
  irtpp::convert_vector(Rindividual_weights, individual_weights);
  irtpp::convert_vector(Rclusters, clusters);
  irtpp::convert_matrix(Rinitial_values, initial_values);

  //Estimation object  
  irtpp::dichotomous::estimation e(Y, dim, model, EMepsilon, 
                                   theta, weights, individual_weights,
                                   clusters, initial_values);

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

  return List::create(Rcpp::Named("zetas") = zetas,
                      Rcpp::Named("Loglikelihood") = e.log_likelihood());;
}


NumericMatrix ltraitscpp ( IntegerMatrix Rdata, unsigned int dim, int model, 
                           NumericMatrix Rzetas,   
                           NumericMatrix Rtheta, NumericVector Rweights, 
                           std::string method,
                           bool by_individuals,
                           NumericMatrix Rinit_traits ) {
  // Converting data types
  irtpp::matrix<char> Y;
  irtpp::matrix<double> zetas;
  irtpp::matrix<double> theta;
  std::vector<double> weights;
  irtpp::matrix<double> init_traits;

  irtpp::convert_matrix(Rdata, Y);
  irtpp::convert_matrix(Rzetas, zetas);
  irtpp::convert_matrix(Rtheta, theta);
  irtpp::convert_vector(Rweights, weights);
  irtpp::convert_matrix(Rinit_traits, init_traits);

  //Estimation object  
  irtpp::dichotomous::estimation e(Y, dim, model, 1e-4, 
                                   theta, weights );
  e.load_multi_initial_values(zetas);

  if ( method == "EAP" ) e.EAP(by_individuals);
  else                   e.MAP(by_individuals);

  //NumericMatrix latent_traits;
  //irtpp::convert_matrix(e.data.latent_traits, latent_traits);

  //Nothing yet
  //TODO, convert from dlib vectors to NumericMatrix
  return NumericMatrix(0, 0, 0);
}