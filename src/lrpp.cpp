#include "lrpp.h"

using namespace Rcpp;

List dichotomous ( IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon,
                   NumericMatrix Rtheta, NumericVector Rweights, 
                   IntegerVector Rindividual_weights, IntegerVector Rclusters,
                   NumericMatrix Rinitial_values ) {
  // Converting data types
  lrpp::matrix<char> Y;
  lrpp::matrix<double> theta;
  std::vector<double> weights;
  std::vector<int> individual_weights;
  std::vector<int> clusters;
  lrpp::matrix<double> initial_values;

  lrpp::convert_matrix(Rdata, Y);
  lrpp::convert_matrix(Rtheta, theta);
  lrpp::convert_vector(Rweights, weights);
  lrpp::convert_vector(Rindividual_weights, individual_weights);
  lrpp::convert_vector(Rclusters, clusters);
  lrpp::convert_matrix(Rinitial_values, initial_values);

  //Estimation object  
  lrpp::dichotomous::estimation e(Y, dim, model, EMepsilon, 
                                   theta, weights, individual_weights,
                                   clusters, initial_values);

  //EM
  e.EMAlgorithm();
  

  NumericMatrix zetas(e.data.p, e.data.d + 2);
  int current_zeta = e.get_iterations() % lrpp::ACCELERATION_PERIOD;

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

List ltraitscpp ( IntegerMatrix Rdata, unsigned int dim, int model, 
                           NumericMatrix Rzetas,   
                           NumericMatrix Rtheta, NumericVector Rweights, 
                           std::string method,
                           bool by_individuals,
                           NumericMatrix Rinit_traits ) {
  // Converting data types
  lrpp::matrix<char> Y;
  lrpp::matrix<double> zetas;
  lrpp::matrix<double> theta;
  std::vector<double> weights;
  lrpp::matrix<double> init_traits;

  lrpp::convert_matrix(Rdata, Y);
  lrpp::convert_matrix(Rzetas, zetas);
  lrpp::convert_matrix(Rtheta, theta);
  lrpp::convert_vector(Rweights, weights);
  lrpp::convert_matrix(Rinit_traits, init_traits);

  

  //Estimation object  
  lrpp::dichotomous::estimation e( Y, dim, model, 1e-4, 
                                   theta, weights );
  e.load_multi_initial_values(zetas);

  if ( method == "EAP" ) e.EAP(by_individuals);
  else                   e.MAP(by_individuals);

  NumericMatrix traits(4, 4, 4);
  return List::create(Rcpp::Named("latent_traits") = traits);

  
  //lrpp::convert_matrix(e.data.latent_traits, traits);

}