#include <Rcpp.h>
#include <vector>
#include "matrix.h"

namespace irtpp {

template <class T>
void convert_matrix ( Rcpp::IntegerMatrix mat, irtpp::matrix<T> &Y ) {
  Y = irtpp::matrix<T>(mat.nrow(), mat.ncol());
  for ( int i = 0; i < mat.nrow(); ++i )
    for ( int j = 0; j < mat.ncol(); ++j )
      Y(i, j) = (T)mat(i,j);
}

template <class T>
void convert_matrix ( Rcpp::NumericMatrix mat, irtpp::matrix<T> &Y ) {
  Y = irtpp::matrix<T>(mat.nrow(), mat.ncol());
  for ( int i = 0; i < mat.nrow(); ++i )
    for ( int j = 0; j < mat.ncol(); ++j )
      Y(i, j) = (T)mat(i,j);
}


template <class T>
void convert_vector ( Rcpp::NumericVector rvector, std::vector<T> &cvector ) {
  cvector = std::vector<T>(rvector.size());
  for ( int i = 0; i < cvector.size(); ++i )
  	cvector[i] = rvector[i];
}

}