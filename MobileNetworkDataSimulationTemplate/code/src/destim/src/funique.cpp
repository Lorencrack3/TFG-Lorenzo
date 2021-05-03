#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// This is a very simple function to remove duplicates in a matrix
// of float numbers. It only scales well when all but a few rows are
// duplicates

// [[Rcpp::export]]
NumericMatrix cppfunique(NumericMatrix mat, double tol) {
  int i, j, k, l;

  j = 1;
  for (i = 1; i < mat.nrow(); ++i) {
    for (k = 0; k < j; ++k) {
      for (l = 0; l < mat.ncol(); ++l)
        if (std::abs(mat(i, l) - mat(k, l)) > tol)
          break;
      if (l == mat.ncol())
        break;
    }
    if (k == j)
      mat(j++, _) = mat(i, _);
  }
  return(mat(Range(0, j - 1), _));
}

// [[Rcpp::export]]
IntegerVector cppfuniqueind(NumericMatrix mat, double tol) {
  int i, j, k, l;
  IntegerVector output(mat.nrow());

  j = 1;
  output(0) = 1;
  for (i = 1; i < mat.nrow(); ++i) {
    for (k = 0; k < j; ++k) {
      for (l = 0; l < mat.ncol(); ++l)
        if (std::abs(mat(i, l) - mat(k, l)) > tol)
          break;
        if (l == mat.ncol())
          break;
    }
    if (k == j) {
      mat(j++, _) = mat(i, _);
      output(k) = i + 1;
    }
  }
  return(output[Range(0, j - 1)]);
}
