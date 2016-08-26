#include "latentregpp.h"

using namespace Rcpp;

List latentregcpp ( IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon,
                   NumericMatrix Rtheta, NumericVector Rweights, 
                   IntegerVector Rindividual_weights, 
                   bool dichotomous_data,
                   IntegerVector Rclusters,
                   NumericMatrix Rinitial_values ) {
  // Converting data types
  latentregpp::matrix<char> Y;
  latentregpp::matrix<double> theta;
  std::vector<double> weights;
  std::vector<int> individual_weights;
  std::vector<int> clusters;
  latentregpp::matrix<double> initial_values;

  latentregpp::convert_matrix(Rdata, Y);
  latentregpp::convert_matrix(Rtheta, theta);
  latentregpp::convert_vector(Rweights, weights);
  latentregpp::convert_vector(Rindividual_weights, individual_weights);
  latentregpp::convert_vector(Rclusters, clusters);
  latentregpp::convert_matrix(Rinitial_values, initial_values);

  if ( dichotomous_data ) {
    //Estimation object  
    latentregpp::dichotomous::estimation e(Y, dim, model, EMepsilon, 
                                     theta, weights, individual_weights,
                                     clusters, initial_values);
    //EM
    e.EMAlgorithm();
    
    NumericMatrix zetas(e.data.p, e.data.d + 2);
    int current_zeta = e.get_iterations() % latentregpp::ACCELERATION_PERIOD;
    int parameters = e.data.m.parameters;

    for ( int i = 0; i < e.data.p; ++i ) {
      int j = 0;
      //a's
      if ( parameters == latentregpp::ONEPL )
        for ( ; j < e.data.d; ++j ) zetas(i, j) = latentregpp::ALPHA_WITH_NO_ESTIMATION;
      else
        for ( ; j < e.data.d; ++j ) zetas(i, j) = e.data.zeta[current_zeta][i](j);
      
      //d
      zetas(i, j) = e.data.zeta[current_zeta][i](j);
      ++j;

      //c
      if ( parameters == latentregpp::THREEPL ) {
        double &c = zetas(i, j);
        c = e.data.zeta[current_zeta][i](j);
        c = 1.0 / (1.0 + exp(-c));
      }
      else 
        zetas(i, j) = 0;  
    }

    return List::create(Rcpp::Named("zetas") = zetas,
                        Rcpp::Named("Loglikelihood") = e.log_likelihood());
  }

  //polytomous data

  //Estimation object  
  latentregpp::polytomous::estimation e(Y, dim, model, EMepsilon, 
                                   theta, weights, individual_weights,
                                   clusters, initial_values);
  //EM
  e.EMAlgorithm();
  
  int max_category = *std::max_element(e.data.categories_item.begin(), e.data.categories_item.end());
  NumericMatrix zetas(e.data.p, e.data.d + max_category - 1);
  std::fill(zetas.begin(), zetas.end(), NumericVector::get_na());
  int current_zeta = e.get_iterations() % latentregpp::ACCELERATION_PERIOD;
  int parameters = e.data.m.parameters;

  for ( int i = 0; i < e.data.p; ++i ) {
    int j = 0;
    //a's
    if ( parameters == latentregpp::ONEPL )
      for ( ; j < e.data.d; ++j ) zetas(i, j) = latentregpp::ALPHA_WITH_NO_ESTIMATION;
    else
      for ( ; j < e.data.d; ++j ) zetas(i, j) = e.data.zeta[current_zeta][i](j);

    //d's
    int categories_item_i = e.data.categories_item[i];
    for ( int h = 1; h < categories_item_i; ++h, ++j )
      zetas(i, j) = e.data.zeta[current_zeta][i](j);
  }

  return List::create(Rcpp::Named("zetas") = zetas,
                      Rcpp::Named("Loglikelihood") = e.log_likelihood());
}

List ltraitscpp ( IntegerMatrix Rdata, unsigned int dim, int model, 
                           NumericMatrix Rzetas,   
                           NumericMatrix Rtheta, NumericVector Rweights, 
                           std::string method,
                           bool by_individuals,
                           bool dichotomous_data,
                           NumericMatrix Rinit_traits ) {
  // Converting data types
  latentregpp::matrix<char> Y;
  latentregpp::matrix<double> zetas;
  latentregpp::matrix<double> theta;
  std::vector<double> weights;
  latentregpp::matrix<double> init_traits;

  latentregpp::convert_matrix(Rdata, Y);
  latentregpp::convert_matrix(Rzetas, zetas);
  latentregpp::convert_matrix(Rtheta, theta);
  latentregpp::convert_vector(Rweights, weights);
  latentregpp::convert_matrix(Rinit_traits, init_traits);
  
  if ( dichotomous_data ) {
    //Estimation object 
    latentregpp::dichotomous::estimation e( Y, dim, model, 1e-4, 
                                     theta, weights );
    e.load_multi_initial_values(zetas);

    //Latent traits
    if ( method == "EAP" ) e.EAP(by_individuals);
    else { 
      if ( init_traits.rows() != 0 ) { 
        //TODO load initial traits
      }
      e.MAP(by_individuals);
    }

    NumericMatrix traits;
    latentregpp::convert_matrix(e.data.latent_traits, traits);
    
    if ( by_individuals ) return List::create(Rcpp::Named("latent_traits") = traits);

    NumericMatrix patterns;
    latentregpp::convert_matrix(e.data.Y, patterns);  
    return List::create(Rcpp::Named("latent_traits") = traits, 
                        Rcpp::Named("patterns") = patterns);
  }

  //Estimation object 
  latentregpp::polytomous::estimation e( Y, dim, model, 1e-4, 
                                   theta, weights );
  e.load_multi_initial_values(zetas);

  //Latent traits
  if ( method == "EAP" ) e.EAP(by_individuals);
  else { 
    if ( init_traits.rows() != 0 ) { 
      //TODO load initial traits
    }
    e.MAP(by_individuals);
  }

  NumericMatrix traits;
  latentregpp::convert_matrix(e.data.latent_traits, traits);
  
  if ( by_individuals ) return List::create(Rcpp::Named("latent_traits") = traits);

  NumericMatrix patterns;
  latentregpp::convert_matrix(e.data.Y, patterns);  
  return List::create(Rcpp::Named("latent_traits") = traits, 
                      Rcpp::Named("patterns") = patterns);
}