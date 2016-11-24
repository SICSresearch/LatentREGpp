#include "latentregpp.h"

using namespace Rcpp;

List itemfitcpp ( IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon,
                   NumericMatrix Rtheta, NumericVector Rweights, 
                   NumericVector Rindividual_weights, 
                   bool dichotomous_data,
                   IntegerVector Rpinned_items,
                   NumericMatrix Rinitial_values,
                   bool verbose ) {
  // Converting data types
  latentregpp::matrix<char> Y;
  latentregpp::matrix<double> theta;
  std::vector<double> weights;
  std::vector<double> individual_weights;
  std::vector<int> pinned_items;
  latentregpp::matrix<double> initial_values;

  latentregpp::convert_matrix(Rdata, Y);
  latentregpp::convert_matrix(Rtheta, theta);
  latentregpp::convert_vector(Rweights, weights);
  latentregpp::convert_vector(Rindividual_weights, individual_weights);
  latentregpp::convert_vector(Rpinned_items, pinned_items);
  latentregpp::convert_matrix(Rinitial_values, initial_values);

  if ( dichotomous_data ) {
    //Estimation object  
    latentregpp::dichotomous::estimation e(Y, dim, model, EMepsilon, 
                                     theta, weights, individual_weights,
                                     pinned_items, initial_values);
    //EM
    e.EMAlgorithm(verbose);
    
    NumericMatrix zetas(e.data.p, e.data.d + 2);
    int current_zeta = e.get_iterations() % latentregpp::ACCELERATION_PERIOD;
    int parameters = e.data.m->parameters;

    for ( int i = 0; i < e.data.p; ++i ) {
      int j = 0;
      //a's and d
      if ( parameters == latentregpp::ONE_PARAMETER ) {
        for ( ; j < e.data.d; ++j ) zetas(i, j) = latentregpp::ALPHA_WITH_NO_ESTIMATION;
        zetas(i, j) = e.data.zeta[current_zeta][i](j - e.data.d);
      }
      else {
        for ( ; j < e.data.d; ++j ) zetas(i, j) = e.data.zeta[current_zeta][i](j);
        zetas(i, j) = e.data.zeta[current_zeta][i](j);
      }
      
      ++j;

      //c
      if ( e.data.m->type == latentregpp::model_type::bayesian && parameters == latentregpp::THREE_PARAMETERS) {
          double &c = zetas(i, j);
          c = e.data.initial_values(i,j);
          c = 1.0 / (1.0 + exp(-c));
      } 
      else if ( parameters == latentregpp::THREE_PARAMETERS ) {
        double &c = zetas(i, j);
        c = e.data.zeta[current_zeta][i](j);
        c = 1.0 / (1.0 + exp(-c));
      }
      else 
        zetas(i, j) = 0;  
    }


    NumericMatrix Rr;
    NumericVector Rf;
    latentregpp::convert_matrix(e.data.r, Rr);
    latentregpp::convert_vector(e.data.f, Rf);

    int h = e.data.m->parameters == 1 ? 1 : e.data.m->parameters - 1 + dim;
    if(e.data.m->parameters == 3) h++;

    double aic = -2 * e.log_likelihood() + 2 * h;
    double bic = -2 * e.log_likelihood() + h * std::log(e.data.N);

    return List::create(Rcpp::Named("zetas") = zetas,
                        Rcpp::Named("Loglik") = e.log_likelihood(),
                        Rcpp::Named("iterations") = e.get_iterations(),
                        Rcpp::Named("r") = Rr,
                        Rcpp::Named("f") = Rf,
                        Rcpp::Named("aic") = aic,
                        Rcpp::Named("bic") = bic);
  }

  //polytomous data

  //Estimation object  
  latentregpp::polytomous::estimation e(Y, dim, model, EMepsilon, 
                                   theta, weights, individual_weights,
                                   pinned_items, initial_values);
  //EM
  e.EMAlgorithm(verbose);
  
  int max_category = *std::max_element(e.data.categories_item.begin(), e.data.categories_item.end());
  NumericMatrix zetas(e.data.p, e.data.d + max_category - 1);
  std::fill(zetas.begin(), zetas.end(), NumericVector::get_na());
  int current_zeta = e.get_iterations() % latentregpp::ACCELERATION_PERIOD;
  int parameters = e.data.m->parameters;

  for ( int i = 0; i < e.data.p; ++i ) {
    int j = 0;
    //a' and d's
    if ( parameters == latentregpp::ONE_PARAMETER ) {
      for ( ; j < e.data.d; ++j ) zetas(i, j) = latentregpp::ALPHA_WITH_NO_ESTIMATION;
      int categories_item_i = e.data.categories_item[i];
      for ( int h = 1; h < categories_item_i; ++h, ++j )
        zetas(i, j) = e.data.zeta[current_zeta][i](j - e.data.d);
    }
    else {
      for ( ; j < e.data.d; ++j ) zetas(i, j) = e.data.zeta[current_zeta][i](j);
      int categories_item_i = e.data.categories_item[i];
      for ( int h = 1; h < categories_item_i; ++h, ++j )
        zetas(i, j) = e.data.zeta[current_zeta][i](j);
    }

    
  }

  int G = e.data.r.size();
  List Rr(G);
  for ( int g = 0; g < G; ++g ) {
    NumericMatrix Rrg;
    latentregpp::convert_matrix(e.data.r[g],  Rrg, max_category);
    Rr[g] = Rrg;
  }

  int h = e.data.m->parameters == 1 ? 1 : e.data.m->parameters - 1 + dim;
  if(e.data.m->parameters == 3) h++;

  double aic = -2 * e.log_likelihood() + 2 * h;
  double bic = -2 * e.log_likelihood() + h * std::log(e.data.N);

  return List::create(Rcpp::Named("zetas") = zetas,
                      Rcpp::Named("Loglik") = e.log_likelihood(),
                      Rcpp::Named("iterations") = e.get_iterations(),
                      Rcpp::Named("r") = Rr,
                      Rcpp::Named("aic") = aic,
                      Rcpp::Named("bic") = bic);
}

List itemfitcpp_bayesian ( IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon,
                   NumericMatrix Rtheta, NumericVector Rweights, 
                   IntegerVector Rindividual_weights, 
                   bool dichotomous_data,
                   IntegerVector Rclusters,
                   NumericMatrix Rinitial_values,
                   bool noguessing,
                   bool verbose ) {
  // Converting data types
  latentregpp::matrix<char> Y;
  latentregpp::matrix<double> theta;
  std::vector<double> weights;
  std::vector<double> individual_weights;
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
                                     clusters, initial_values, true, noguessing);
    //EM
    e.EMAlgorithm(verbose);
    
    NumericMatrix zetas(e.data.p, e.data.d + 2);
    int current_zeta = e.get_iterations() % latentregpp::ACCELERATION_PERIOD;
    int parameters = e.data.m->parameters;

    for ( int i = 0; i < e.data.p; ++i ) {
      int j = 0;
      //a's and d
      if ( parameters == latentregpp::ONE_PARAMETER ) {
        for ( ; j < e.data.d; ++j ) zetas(i, j) = latentregpp::ALPHA_WITH_NO_ESTIMATION;
        zetas(i, j) = e.data.zeta[current_zeta][i](j - e.data.d);
      }
      else {
        for ( ; j < e.data.d; ++j ) zetas(i, j) = e.data.zeta[current_zeta][i](j);
        zetas(i, j) = e.data.zeta[current_zeta][i](j);
      }
      
      ++j;

      //c
      if ( e.data.m->type == latentregpp::model_type::bayesian && parameters == latentregpp::THREE_PARAMETERS && e.data.noguessing) {
          double &c = zetas(i, j);
          c = e.data.initial_values(i,j);
          c = 1.0 / (1.0 + exp(-c));
      } 
      else if ( parameters == latentregpp::THREE_PARAMETERS ) {
        double &c = zetas(i, j);
        c = e.data.zeta[current_zeta][i](j);
        c = 1.0 / (1.0 + exp(-c));
      }
      else 
        zetas(i, j) = 0;  
    }


    NumericMatrix Rr;
    NumericVector Rf;
    latentregpp::convert_matrix(e.data.r, Rr);
    latentregpp::convert_vector(e.data.f, Rf);

    return List::create(Rcpp::Named("zetas") = zetas,
                        Rcpp::Named("Loglik") = e.log_likelihood(),
                        Rcpp::Named("iterations") = e.get_iterations(),
                        Rcpp::Named("r") = Rr,
                        Rcpp::Named("f") = Rf);
  }

}

List personfitcpp ( IntegerMatrix Rdata, unsigned int dim, int model, 
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

  latentregpp::convert_matrix(Rdata, Y);
  latentregpp::convert_matrix(Rzetas, zetas);
  latentregpp::convert_matrix(Rtheta, theta);
  latentregpp::convert_vector(Rweights, weights);
  
  if ( dichotomous_data ) {
    //Estimation object
    latentregpp::dichotomous::estimation e( Y, dim, model, 1e-4, 
                                     theta, weights );
    e.load_multi_initial_values(zetas);

    //Latent traits
    if ( method == "EAP" ) e.EAP(by_individuals);
    else { 
      latentregpp::convert_matrix(Rinit_traits, e.data.latent_traits);
      e.MAP(by_individuals);
    }

    NumericMatrix traits;
    latentregpp::convert_matrix(e.data.latent_traits, traits);
    
    if ( by_individuals ) return List::create(Rcpp::Named("latent_traits") = traits);

    NumericMatrix patterns;
    IntegerMatrix freq;
    latentregpp::convert_matrix(e.data.Y, patterns);  
    latentregpp::convert_vector(e.data.nl, freq);
    return List::create(Rcpp::Named("latent_traits") = traits, 
                        Rcpp::Named("patterns") = patterns,
                        Rcpp::Named("freqs") = freq);
  }

  //Estimation object 
  //printf("here 1");
  latentregpp::polytomous::estimation e( Y, dim, model, 1e-4, 
                                   theta, weights );
  //printf("here 2");
  e.load_multi_initial_values(zetas); //Revisar TODO
  //printf("here 3");
  //Latent traits
  if ( method == "EAP" ) e.EAP(by_individuals);
  else { 
    latentregpp::convert_matrix(Rinit_traits, e.data.latent_traits);
    e.MAP(by_individuals);
  }

  NumericMatrix traits;
  latentregpp::convert_matrix(e.data.latent_traits, traits);
  
  if ( by_individuals ) return List::create(Rcpp::Named("latent_traits") = traits);

  NumericMatrix patterns;
  IntegerMatrix freq;
  latentregpp::convert_matrix(e.data.Y, patterns);  
  latentregpp::convert_vector(e.data.nl, freq);

  return List::create(Rcpp::Named("latent_traits") = traits, 
                      Rcpp::Named("patterns") = patterns,
                      Rcpp::Named("freqs") = freq);
}