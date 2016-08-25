#include <Rcpp.h>
#include <vector>
#include "matrix.h"
#include "constants.h"

namespace lrpp {

/**
* Helper functions to convert R data to C++ data and viceversa
*/

/**
 * Converts from Rcpp::IntegerMatrix to lrpp::matrix<T>
 * @param mat1 source object
 * @param mat2 destination object
 */
template <class T>
void convert_matrix ( Rcpp::IntegerMatrix &mat1, lrpp::matrix<T> &mat2 ) {
  mat2 = lrpp::matrix<T>(mat1.nrow(), mat1.ncol());
  for ( int i = 0; i < mat1.nrow(); ++i )
    for ( int j = 0; j < mat1.ncol(); ++j )
      mat2(i, j) = (T)mat1(i,j);
}

/**
 * Converts from Rcpp::NumericMatrix to lrpp::matrix<T>
 * @param mat1 source object
 * @param mat2 destination object
 */
template <class T>
void convert_matrix ( Rcpp::NumericMatrix &mat1, lrpp::matrix<T> &mat2 ) {
  mat2 = lrpp::matrix<T>(mat1.nrow(), mat1.ncol());
  for ( int i = 0; i < mat1.nrow(); ++i )
    for ( int j = 0; j < mat1.ncol(); ++j )
      mat2(i, j) = (T)mat1(i,j);
}

/**
 * Converts from lrpp::matrix<T> to Rcpp::NumericMatrix
 * @param mat1 source object
 * @param mat2 destination object
 */
template <class T>
void convert_matrix ( lrpp::matrix<T> &mat1, Rcpp::NumericMatrix &mat2 ) {
  mat2 = Rcpp::NumericMatrix(mat1.rows(), mat1.columns(0));
  for ( size_t i = 0; i < mat1.rows(); ++i ) 
    for ( int j = 0; j < mat1.columns(i); ++j )
      mat2(i, j) = mat1(i, j);
}

/**
 * Converts from std::vector<lrpp::optimizer_vector> to Rcpp::NumericMatrix
 * @param mat1 source object
 * @param mat2 destination object
 */
void convert_matrix ( std::vector<lrpp::optimizer_vector> &mat1, Rcpp::NumericMatrix &mat2 ) {
  mat2 = Rcpp::NumericMatrix(int(mat1.size()), int(mat1[0].size()));
  for ( size_t i = 0; i < mat1.size(); ++i ) 
    for ( int j = 0; j < mat1[i].size(); ++j )
      mat2(i, j) = mat1[i](j);
}

/**
 * Converts from Rcpp::NumericVector to std::vector<T> 
 * @param rvector source object
 * @param cvector destination object
 */
template <class T>
void convert_vector ( Rcpp::NumericVector &rvector, std::vector<T> &cvector ) {
  cvector = std::vector<T>(rvector.size());
  for ( int i = 0; i < cvector.size(); ++i )
    cvector[i] = rvector[i];
}

/**
 * Converts from Rcpp::IntegerVector to std::vector<T> 
 * @param rvector source object
 * @param cvector destination object
 */
template <class T>
void convert_vector ( Rcpp::IntegerVector &rvector, std::vector<T> &cvector ) {
  cvector = std::vector<T>(rvector.size());
  for ( int i = 0; i < cvector.size(); ++i )
    cvector[i] = rvector[i];
}

}