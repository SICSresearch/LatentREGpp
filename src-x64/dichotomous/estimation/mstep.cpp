/*
 * mstep.cpp
 *
 *  Created on: 13/04/2016
 *      Author: Milder
 */

#include "../../dichotomous/estimation/mstep.h"

namespace latentregpp {

namespace dichotomous {

Qi::Qi (int i, estimation_data *d, int current_zeta_obj) : i(i), data(d), current_zeta_obj(current_zeta_obj) { }

double Qi::operator() ( const optimizer_vector& item_i ) const {
  //Value of Qi
  double value = 0;
  //Number of quadrature points
  int G = data->G;
  //Matrix r
  matrix<double> &r = data->r;
  //Latent trait vectors
  matrix<double> &theta = data->theta;
  //f
  std::vector<double> &f = data->f;
  
  
  //For Bayessian
  std::vector<optimizer_vector> &current_zeta = data->zeta[current_zeta_obj];
  int Nind = data->N;
  
  
  for ( int g = 0; g < G; ++g ) {
    std::vector<double> &theta_g = *theta.get_pointer_row(g);
    double P_gi = data->m->P(theta_g, item_i);
    value += r(g, i) * log(P_gi) + (f[g] - r(g, i)) * log(1 - P_gi);
  }
  
  //Log(Pzetai)
  //Bayessian mode
  if(false) {
    
    double coef = -(Nind/2);
    
    double pzetai = 0;
    bool guessing_parameter = data->m->parameters == THREE_PARAMETERS;
    
    for (int j = 0; j < current_zeta[i].size()-guessing_parameter; ++j)
    {
      double miu_alpha = 1;
      double sigma_alpha = 1;
      
      if(j+1==current_zeta[i].size()-guessing_parameter) {
        double miu_d = 0;
        double sigma_d = 4;
        pzetai += pow((current_zeta[i](j)-miu_d),2) / pow(sigma_d,2);
        break;			
      }
      
      pzetai += pow((current_zeta[i](j)-miu_alpha),2) / pow(sigma_alpha,2);
    }
    
    if(guessing_parameter) {
      double miu_gamma = 0.01;//-4.59512;
      double sigma_gamma = 0.0009;//7;	 
      pzetai += pow((current_zeta[i](current_zeta[i].size() - 1)-miu_gamma),2) / pow(sigma_gamma,2);	
    }
    
    pzetai *= coef;
    pzetai += value;
    return pzetai;
  }
  
  return value;
}

double Mstep(estimation_data &data, int current) {
  double max_difference = 0.0;
  int next = (current + 1) % ACCELERATION_PERIOD;
  
  int &p = data.p;
  std::vector<optimizer_vector> &current_zeta = data.zeta[current];
  std::vector<optimizer_vector> &next_zeta = data.zeta[next];
  
  if ( next_zeta.empty() ) {
    for ( auto c : current_zeta )
      next_zeta.push_back(c);
  }
  
  std::set<int> &pinned_items = data.pinned_items;
  
  /**
  * Log likelihood must be optimized for every item
  * */
  
  #pragma omp parallel for schedule(dynamic) reduction(max:max_difference)
  for ( int i = 0; i < p; ++i ) {
    /**
    * If it is multidimensional and this is one of the pinned items
    * i.e the first item of a dimension
    * this item is just skipped
    * */
    
    if ( pinned_items.count(i) ) continue;
    
    
    next_zeta[i] = current_zeta[i];
    
    //Calling BFGS from dlib to optimize Qi with approximate derivatives
    dlib::find_max_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                                                 dlib::objective_delta_stop_strategy(OPTIMIZER_DELTA_STOP),
                                                 Qi(i, &data, current), next_zeta[i], -1);

    //for ( int j = 0; j < next_zeta[i].size(); ++j )
      //Rprintf("%lf %lf\n", current_zeta[i](j), next_zeta[i](j));
    
    //Computing difference of current item
    if ( data.m->parameters < THREE_PARAMETERS ) {
      for ( int j = 0; j < next_zeta[i].size(); ++j )
        max_difference = std::max(max_difference, std::abs(next_zeta[i](j) - current_zeta[i](j)));
    } else {
      for ( int j = 0; j < next_zeta[i].size() - 1; ++j )
        max_difference = std::max(max_difference, std::abs(next_zeta[i](j) - current_zeta[i](j)));
      double c_current = current_zeta[i](current_zeta[i].size() - 1);
      double c_next = next_zeta[i](next_zeta[i].size() - 1);
      
      c_current = 1.0 / (1.0 + exp(-c_current));
      c_next = 1.0 / (1.0 + exp(-c_next));
      
      max_difference = std::max(max_difference, std::abs(c_next - c_current));
    }
  }
  
  return max_difference;
}

}

} /* namespace latentregpp */
