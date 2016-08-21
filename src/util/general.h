#include <Rcpp.h>
#include <vector>
#include "matrix.h"

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
