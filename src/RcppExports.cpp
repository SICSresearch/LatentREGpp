// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dichotomous
List dichotomous(IntegerMatrix Rdata, unsigned int dim, int model, double EMepsilon, NumericMatrix theta, NumericVector weights, IntegerVector Rclusters, std::string initial_values);
RcppExport SEXP MulTRI_dichotomous(SEXP RdataSEXP, SEXP dimSEXP, SEXP modelSEXP, SEXP EMepsilonSEXP, SEXP thetaSEXP, SEXP weightsSEXP, SEXP RclustersSEXP, SEXP initial_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Rdata(RdataSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< double >::type EMepsilon(EMepsilonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Rclusters(RclustersSEXP);
    Rcpp::traits::input_parameter< std::string >::type initial_values(initial_valuesSEXP);
    __result = Rcpp::wrap(dichotomous(Rdata, dim, model, EMepsilon, theta, weights, Rclusters, initial_values));
    return __result;
END_RCPP
}
