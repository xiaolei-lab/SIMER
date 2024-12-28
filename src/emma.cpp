#include "simer.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]

template <typename T>
arma::vec BigRowMean(XPtr<BigMatrix> pMat, int threads = 0){

  omp_setup(threads);

  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

  int ind = pMat->ncol();
  int j, k, m = pMat->nrow();
  double p1 = 0.0;
  arma::vec mean(m);

  #pragma omp parallel for private(p1, k)
  for (j = 0; j < m; j++) {
    p1 = 0.0;
    for (k = 0; k < ind; k++) {
      p1 += bigm[k][j];
    }
    mean[j] = p1 / ind;
  }

  return mean;
}

arma::vec BigRowMean(SEXP pBigMat, int threads = 0){
  
  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()) {
  case 1:
    return BigRowMean<char>(xpMat, threads);
  case 2:
    return BigRowMean<short>(xpMat, threads);
  case 4:
    return BigRowMean<int>(xpMat, threads);
  case 8:
    return BigRowMean<double>(xpMat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template <typename T>
arma::mat emma_kinship(XPtr<BigMatrix> pMat, int threads = 0, bool verbose=true) {
  omp_setup(threads);
  
  int i, j, k, m = pMat->nrow(), n = pMat->ncol();
  double s;

  arma::mat K(n, n, fill::ones);

  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

  arma::vec rowMeans = BigRowMean(pMat, threads);

  MinimalProgressBar pb;
  Progress p(n, verbose, pb);

  if (verbose) { Rcout << " Computing EMMA Kinship Matrix..." << endl; }

  #pragma omp parallel for private(i, j, k, s)
  for (i = 0; i < n; i++) {
    for (j = i+1; j < n; j++) {
      s = 0;
      for (k = 0; k < m; k++) {
        if (bigm[i][k] == bigm[j][k]) {
          s = s + 1;
        } else {
          if ((bigm[i][k] == 1) || (bigm[j][k] == 1)) {
            if (rowMeans[k] == 1) {
              if ((bigm[i][k] + bigm[j][k]) == 1) {
                s = s + 1;
              }
            } else {
              s = s + 0.5;
            }
          }
        }
      }
      K(i, j) = s / m;
      K(j, i) = K(i, j);
    }
    if ( ! Progress::check_abort() ) { p.increment(); }
  }

  return K;
}

// [[Rcpp::export]]
arma::mat emma_kinship(SEXP pBigMat, int threads = 0, bool verbose=true){
  
  XPtr<BigMatrix> xpMat(pBigMat);

  switch(xpMat->matrix_type()) {
  case 1:
    return emma_kinship<char>(xpMat, threads, verbose);
  case 2:
    return emma_kinship<short>(xpMat, threads, verbose);
  case 4:
    return emma_kinship<int>(xpMat, threads, verbose);
  case 8:
    return emma_kinship<double>(xpMat, threads, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}
