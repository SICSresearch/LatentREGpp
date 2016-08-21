#include <Rcpp.h>
#include "simulation/simulation.h"
using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::plugins(openmp)]]
template <class T>
void convert_matrix ( Rcpp::IntegerMatrix mat, irtpp::matrix<T> &Y ) {
  Y = irtpp::matrix<T>(mat.nrow(), mat.ncol());
  for ( int i = 0; i < mat.nrow(); ++i )
    for ( int j = 0; j < mat.ncol(); ++j )
      Y(i, j) = (T)mat(i,j);
}

void convert_vector(Rcpp::NumericVector input, std::vector<int> &output){
  output = std::vector<int>(input.size());
  for(int i = 0; i < input.size(); i++)
    output.push_back(output[i]);
}

//' Estimates the parameters of a test
//'
//' @param RData Input dataset.
//' @export
// [[Rcpp::export]]
List multiTest_dico(IntegerMatrix RData,unsigned int dim,std::string wd){
  irtpp::matrix<char> Y;
  convert_matrix(RData, Y);
  Rcout << "Data Loaded" << std::endl;
  Rcout << "Data Rows: " << Y.rows() <<std::endl;
  Rcout << "Data Columns: " << Y.get_row(1).size() <<std::endl;
  Rcout << "Dimension: " << dim << std::endl;

  //std::vector<int> group_items;
  //convert_vector(no_items,group_items);
  irtpp::setwd(wd);
  Rcout << "Working Directory: " << irtpp::getwd() << "\nInput: "<< wd << std::endl;
  irtpp::dichotomous::estimation e(Y,dim);

  //convert_matrix(quad, e.estimation_data.theta);
  // convert_matrix(weights, e.estimation_data.weights);
  // e.estimation_data.G = G;
  // e.build_matrixes();

  Rcout << "Estimation Object Created" << std::endl;
  e.EMAlgorithm();
  Rcout << "Estimation Completed" << std::endl;

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
