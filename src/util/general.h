#include <Rcpp.h>
#include <vector>
#include "matrix.h"
#include "constants.h"

namespace lrpp {

template <class T>
void convert_matrix ( Rcpp::IntegerMatrix &mat, lrpp::matrix<T> &Y ) {
  Y = lrpp::matrix<T>(mat.nrow(), mat.ncol());
  for ( int i = 0; i < mat.nrow(); ++i )
    for ( int j = 0; j < mat.ncol(); ++j )
      Y(i, j) = (T)mat(i,j);
}

template <class T>
void convert_matrix ( Rcpp::NumericMatrix &mat, lrpp::matrix<T> &Y ) {
  Y = lrpp::matrix<T>(mat.nrow(), mat.ncol());
  for ( int i = 0; i < mat.nrow(); ++i )
    for ( int j = 0; j < mat.ncol(); ++j )
      Y(i, j) = (T)mat(i,j);
}

void convert_matrix ( std::vector<lrpp::optimizer_vector> &mat1, Rcpp::NumericMatrix &mat2 ) {
  mat2 = Rcpp::NumericMatrix(int(mat1.size()), int(mat1[0].size()));
  for ( size_t i = 0; i < mat1.size(); ++i ) 
    for ( int j = 0; j < mat1[i].size(); ++j )
      mat2(i, j) = mat1[i](j);
}

template <class T>
void convert_vector ( Rcpp::NumericVector &rvector, std::vector<T> &cvector ) {
  cvector = std::vector<T>(rvector.size());
  for ( int i = 0; i < cvector.size(); ++i )
    cvector[i] = rvector[i];
}

template <class T>
void convert_vector ( Rcpp::IntegerVector &rvector, std::vector<T> &cvector ) {
  cvector = std::vector<T>(rvector.size());
  for ( int i = 0; i < cvector.size(); ++i )
    cvector[i] = rvector[i];
}

}