#include <RcppArmadillo.h>
#include "simer_omp.h"
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <progress.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

template<typename T>
NumericVector FilterMind(XPtr<BigMatrix> pMat, double NA_C, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  
  size_t i, j, m = pMat->nrow(), n = pMat->ncol();
  NumericVector colNumNA(n, 0);
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++) {
      if (bigm[j][i] == NA_C) { colNumNA[j] += 1;  }
    }
  }
  
  return colNumNA;
}

NumericVector FilterMind(const SEXP pBigMat, double NA_C, int threads=0) {
    XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return FilterMind<char>(xpMat, NA_CHAR, threads);
  case 2:
    return FilterMind<short>(xpMat, NA_SHORT, threads);
  case 4:
    return FilterMind<int>(xpMat, NA_INTEGER, threads);
  case 8:
    return FilterMind<double>(xpMat, NA_REAL, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
NumericVector FilterGeno(XPtr<BigMatrix> pMat, double NA_C, IntegerVector rowIdx, IntegerVector colIdx, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  size_t i, j;
  NumericVector rowNumNA(rowIdx.size(), 0);
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < colIdx.size(); j++) {
    for (i = 0; i < rowIdx.size(); i++) {
      if (bigm[colIdx[j]][rowIdx[i]] == NA_C) { rowNumNA[i] += 1;  }
    }
  }
  
  return rowNumNA;
}

NumericVector FilterGeno(const SEXP pBigMat, double NA_C, IntegerVector rowIdx, IntegerVector colIdx, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return FilterGeno<char>(xpMat, NA_CHAR, rowIdx, colIdx, threads);
  case 2:
    return FilterGeno<short>(xpMat, NA_SHORT, rowIdx, colIdx, threads);
  case 4:
    return FilterGeno<int>(xpMat, NA_INTEGER, rowIdx, colIdx, threads);
  case 8:
    return FilterGeno<double>(xpMat, NA_REAL, rowIdx, colIdx, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
arma::mat CalGenoFreq(XPtr<BigMatrix> pMat, double NA_C, IntegerVector rowIdx, IntegerVector colIdx, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  
  size_t i, j;
  arma::mat genoFreq(rowIdx.size(), 3, fill::zeros);

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (i = 0; i < rowIdx.size(); i++) {
    for (j = 0; j < colIdx.size(); j++) {
      if (bigm[colIdx[j]][rowIdx[i]] == 0) {
        genoFreq(i, 0) = genoFreq(i, 0) + 1; 
      } else if (bigm[colIdx[j]][rowIdx[i]] == 1) {
        genoFreq(i, 1) = genoFreq(i, 1) + 1;
      } else if (bigm[colIdx[j]][rowIdx[i]] == 2) {
        genoFreq(i, 2) = genoFreq(i, 2) + 1;
      }
    }
  }
  
  return genoFreq;
}

arma::mat CalGenoFreq(const SEXP pBigMat, IntegerVector rowIdx, IntegerVector colIdx, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return CalGenoFreq<char>(xpMat, NA_CHAR, rowIdx, colIdx, threads);
  case 2:
    return CalGenoFreq<short>(xpMat, NA_SHORT, rowIdx, colIdx, threads);
  case 4:
    return CalGenoFreq<int>(xpMat, NA_INTEGER, rowIdx, colIdx, threads);
  case 8:
    return CalGenoFreq<double>(xpMat, NA_REAL, rowIdx, colIdx, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

NumericVector FilterMAF(arma::mat genoFreq, int threads=0) {
  omp_setup(threads);
  
  IntegerVector freq0 = wrap(genoFreq.col(0));
  IntegerVector freq1 = wrap(genoFreq.col(1));
  IntegerVector freq2 = wrap(genoFreq.col(2));

  size_t i;
  NumericVector MAF(genoFreq.n_rows); MAF.fill(0);
  
  #pragma omp parallel for schedule(dynamic) private(i)
  for (i = 0; i < genoFreq.n_rows; i++) {
    MAF[i] = (freq0[i] + freq1[i] * 0.5) / 
      (freq0[i]+ freq1[i] + freq2[i]);
    MAF[i] = MAF[i] <= 0.5 ? MAF[i] : (1 - MAF[i]);
  }
  
  return MAF;
}

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) {
  if ((obs_hom1 < 0) || (obs_hom2 < 0) || (obs_hets < 0)) {
    Rcpp::stop("FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a negative count", obs_hets, obs_hom1, obs_hom2);
  }
  
  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
  
  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes   = obs_hets + obs_homc + obs_homr;
  
  double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
  if (het_probs == NULL) {
    Rcpp::stop("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities");
  }
  
  int i;
  for (i = 0; i <= rare_copies; i++) {
    het_probs[i] = 0.0;
  }
  
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
  
  if ((rare_copies & 1) ^ (mid & 1)) {
    mid++;
  }
  
  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;
  
  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
    / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];
    
    curr_homr++;
    curr_homc++;
  }
  
  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
    /((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];
    
    curr_homr--;
    curr_homc--;
  }
  
  for (i = 0; i <= rare_copies; i++) {
    het_probs[i] /= sum;
  }
  
  double p_hwe = 0.0;
  
  for (i = 0; i <= rare_copies; i++) {
    if (het_probs[i] > het_probs[obs_hets]) {
      continue;
    }
    p_hwe += het_probs[i];
  }
  
  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
  
  free(het_probs);
  
  return p_hwe;
}

NumericVector FilterHWE(arma::mat genoFreq, int threads=0) {
  omp_setup(threads);
  
  size_t i;
  IntegerVector freq0 = wrap(genoFreq.col(0));
  IntegerVector freq1 = wrap(genoFreq.col(1));
  IntegerVector freq2 = wrap(genoFreq.col(2));
  NumericVector PVAL(genoFreq.n_rows); PVAL.fill(0);
  
  #pragma omp parallel for schedule(dynamic) private(i)
  for (i = 0; i < genoFreq.n_rows; i++) {
    PVAL[i] = SNPHWE(freq1[i], freq0[i], freq2[i]);
  }
  
  return PVAL;
}

template<typename T>
List GenoFilter(XPtr<BigMatrix> pMat, double NA_C, Nullable<IntegerVector> keepInds=R_NilValue, Nullable<double> filterGeno=R_NilValue, Nullable<double> filterHWE=R_NilValue, Nullable<double> filterMind=R_NilValue, Nullable<double> filterMAF=R_NilValue, int threads=0, bool verbose=true) {

  double m = pMat->nrow(), n = pMat->ncol();
  IntegerVector  keepRows = seq(0, m - 1);
  IntegerVector keepCols;
  if (keepInds.isNull()) {
    keepCols = seq(0, n - 1);
  } else {
    keepCols = as<IntegerVector>(keepInds);
    keepCols = keepCols - 1;
    n = keepCols.size();
  }
  
  double fgeno = 0, fhwe = 0, fmaf = 0, fmind = 0;
  if (filterGeno.isNotNull()) { fgeno = as<double>(filterGeno); }
  if (filterHWE.isNotNull() ) { fhwe  = as<double>(filterHWE ); }
  if (filterMAF.isNotNull() ) { fmaf  = as<double>(filterMAF); }
  if (filterMind.isNotNull()) { fmind = as<double>(filterMind); }
  
  if (verbose) {
    Rcout << " Options in effect:" << endl;
    if (keepInds.isNotNull()  ) { Rcout << "   --keep-ind filePed "       << endl; }
    if (filterGeno.isNotNull()) { Rcout << "   --geno " << fgeno  << endl; }
    if (filterHWE.isNotNull() ) { Rcout << "   --hwe "  << fhwe   << endl; }
    if (filterMAF.isNotNull() ) { Rcout << "   --maf "  << fmaf   << endl; }
    if (filterMind.isNotNull()) { Rcout << "   --mind " << fmind  << endl; }
    Rcout << endl;
    Rcout << " Detect " << n << " samples and " << m << " variants" << endl;
    Rcout << endl;
  }
  
  if (filterMind.isNotNull()) {
    if (verbose) { Rcout << " Calculating sample missingness rates..."; }
    NumericVector colNumNA = FilterMind(pMat, NA_C, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepCols = keepCols[colNumNA/m < fmind];
    if (verbose) {
      Rcout << " " << (n - keepCols.size())  << " samples removed due to missing genotype data (--mind)." << endl;
      n = keepCols.size();
      Rcout << " " << n << " samples remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  if (filterGeno.isNotNull()) {
    if (verbose) { Rcout << " Calculating variant missingness rates..."; }
    NumericVector rowNumNA = FilterGeno(pMat, NA_C, keepRows, keepCols, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepRows = keepRows[rowNumNA/n < fgeno];
    if (verbose) {
      Rcout << " " << (m - keepRows.size()) << " variants removed due to missing genotype data (--geno)." << endl;
      m = keepRows.size();
      Rcout << " " << m << " variants remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  arma::mat genoFreq;
  if (filterMAF.isNotNull() || filterHWE.isNotNull()) {
    if (verbose) { Rcout << " Calculating Genotype Frequencies..."; }
    genoFreq = CalGenoFreq(pMat, keepRows, keepCols, threads);
    if (verbose) {  Rcout << " done." << endl << endl; }
  }
  
  if (filterHWE.isNotNull()) {
    if (verbose) { Rcout << " Performing Hardy-Weinberg test..."; }
    NumericVector PVAL = FilterHWE(genoFreq, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepRows = keepRows[PVAL > fhwe];
    arma::vec armaPVAL = as<arma::vec>(PVAL);
    genoFreq = genoFreq.rows(arma::find(armaPVAL > fhwe));
    if (verbose) {
      Rcout << " " << (m - keepRows.size()) << " variants removed due to exceeding HWE-P-Value (--hwe)." << endl;
      m = keepRows.size();
      Rcout << " " << m << " variants remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  if (filterMAF.isNotNull()) {
    if (verbose) { Rcout << " Calculating Minor Allele Frequencies..."; }
    NumericVector MAF = FilterMAF(genoFreq, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepRows = keepRows[MAF >= fmaf];
    if (verbose) {
      Rcout << " " << (m - keepRows.size()) << " variants removed due to exceeding MAF (--maf)." << endl;
      m = keepRows.size();
      Rcout << " " << m << " variants remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  keepRows = keepRows + 1;
  keepCols = keepCols + 1;
  List genoInfo = List::create(Named("keepRows") = keepRows,
                                   _["keepCols"] = keepCols);
  return genoInfo;
}

// [[Rcpp::export]]
List GenoFilter(const SEXP pBigMat, Nullable<IntegerVector> keepInds=R_NilValue, Nullable<double> filterGeno=R_NilValue, Nullable<double> filterHWE=R_NilValue, Nullable<double> filterMind=R_NilValue, Nullable<double> filterMAF=R_NilValue, int threads=0, bool verbose=true) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return GenoFilter<char>(xpMat, NA_CHAR, keepInds, filterGeno, filterHWE, filterMind, filterMAF, threads, verbose);
  case 2:
    return GenoFilter<short>(xpMat, NA_SHORT, keepInds, filterGeno, filterHWE, filterMind, filterMAF, threads, verbose);
  case 4:
    return GenoFilter<int>(xpMat, NA_INTEGER, keepInds, filterGeno, filterHWE, filterMind, filterMAF, threads, verbose);
  case 8:
    return GenoFilter<double>(xpMat, NA_REAL, keepInds, filterGeno, filterHWE, filterMind, filterMAF, threads, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void Mat2BigMat(XPtr<BigMatrix> pMat, IntegerMatrix mat, Nullable<IntegerVector> colIdx=R_NilValue, int op=1, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  
  IntegerVector ci;
  if (colIdx.isNull()) {
    ci = seq(0, mat.ncol() - 1);
  } else {
    ci = as<IntegerVector>(colIdx);
    ci = ci - 1;
  }
  
  size_t i, j, m = mat.nrow(), n = ci.length();
  op = op - 1;
  if (m != pMat->nrow()) {
    Rcpp::stop("'bigmat' and 'mat' should have the same marker number!");
  }
  if (op + n > pMat->ncol()) {
    Rcpp::stop("'mat' cannot be intert to bigmat completely!");
  }
  if (max(ci) + 1 > mat.ncol()) {
    Rcpp::stop("'colIdx' is out of bound!");
  }
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++) {
      bigmat[op + j][i] = mat(i, ci[j]);
    }
  }

}

// [[Rcpp::export]]
void Mat2BigMat(const SEXP pBigMat, IntegerMatrix mat, Nullable<IntegerVector> colIdx=R_NilValue, int op=1, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return Mat2BigMat<char>(xpMat, mat, colIdx, op, threads);
  case 2:
    return Mat2BigMat<short>(xpMat, mat, colIdx, op, threads);
  case 4:
    return Mat2BigMat<int>(xpMat, mat, colIdx, op, threads);
  case 8:
    return Mat2BigMat<double>(xpMat, mat, colIdx, op, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void BigMat2BigMat(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pmat, Nullable<NumericVector> colIdx=R_NilValue, int op=1, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pmat);
  
  NumericVector ci;
  if (colIdx.isNull()) {
    ci = seq(0, pmat->ncol() - 1);
  } else {
    ci = as<NumericVector>(colIdx);
    ci = ci - 1;
  }
  
  size_t i, j, m = pmat->nrow(), n = ci.length();
  op = op - 1;
  if (m != pMat->nrow()) {
    Rcpp::stop("'bigmat' and 'pmat' should have the same marker number!");
  }
  if (op + n > pMat->ncol()) {
    Rcpp::stop("'pmat' cannot be intert to bigmat completely!");
  }
  if (max(ci) + 1 > pmat->ncol()) {
    Rcpp::stop("'colIdx' is out of bound!");
  }
  
  IntegerMatrix mat(pmat->nrow(), pmat->ncol());
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < pmat->ncol(); j++) {
    for (i = 0; i < m; i++) {
      mat(i, j) = bigm[j][i];
    }
  }
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++) {
      bigmat[op + j][i] = mat(i, ci[j]);
    }
  }
  
}

// [[Rcpp::export]]
void BigMat2BigMat(const SEXP pBigMat, const SEXP pBigmat, Nullable<NumericVector> colIdx=R_NilValue, int op=1, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  XPtr<BigMatrix> xpmat(pBigmat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return BigMat2BigMat<char>(xpMat, xpmat, colIdx, op, threads);
  case 2:
    return BigMat2BigMat<short>(xpMat, xpmat, colIdx, op, threads);
  case 4:
    return BigMat2BigMat<int>(xpMat, xpmat, colIdx, op, threads);
  case 8:
    return BigMat2BigMat<double>(xpMat, xpmat, colIdx, op, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void GenoMixer(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pmat, IntegerVector sirIdx, IntegerVector damIdx, int nBlock=100, int op=1, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pmat);
  
  sirIdx = sirIdx - 1;
  damIdx = damIdx - 1;
  
  size_t op_row, ed_row, i, j, k, m, n, judpar, kidIdx;
  std::random_device rd;
  m = pmat->nrow(); 
  n = damIdx.length();
  op = op - 1;
  
  if (m != pMat->nrow()) {
    Rcpp::stop("'bigmat' and 'pmat' should have the same marker number!");
  }
  if (op + n > pMat->ncol()) {
    Rcpp::stop("'pmat' cannot be intert to bigmat completely!");
  }
  if ((max(sirIdx) > pmat->ncol()) || (max(damIdx) > pmat->ncol())) {
    Rcpp::stop("'sirIdx' or 'damIdx' is out of bound!");
  }
  if (sirIdx.length() != damIdx.length()) {
    Rcpp::stop("'sirIdx' and 'damIdx' should have the same length!");
  }
  
  int len_block = m / nBlock;
  int tail_block = m % nBlock + len_block;
  IntegerVector nInblock(nBlock);
  IntegerVector accum_block(nBlock);
  for (i = 0; i < nBlock; i++) {
    nInblock[i] = len_block;
  }
  nInblock[nBlock - 1] = tail_block;
  accum_block[0] = nInblock[0];
  for (i = 1; i < nBlock; i++) {
    accum_block[i] = accum_block[i - 1] + nInblock[i];
  }
  
  IntegerMatrix mat(pmat->nrow(), pmat->ncol());
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < pmat->ncol(); j++) {
    for (i = 0; i < m; i++) {
      mat(i, j) = bigm[j][i];
    }
  }
  
  #pragma omp parallel for schedule(dynamic) private(i, j, k, op_row, ed_row, judpar, kidIdx)
  for (k = 0; k < nBlock; k++) {
    ed_row = accum_block[k];
    op_row = ed_row - nInblock[k];
    for (j = 0; j < n; j++) {
      judpar = rd();
      kidIdx = judpar % 2 == 0 ? sirIdx[j] : damIdx[j];
      for (i = op_row; i < ed_row; i++) {
        bigmat[op + j][i] = mat(i, kidIdx);
      }
    }
  }
  
}

// [[Rcpp::export]]
void GenoMixer(const SEXP pBigMat, const SEXP pBigmat, IntegerVector sirIdx, IntegerVector damIdx, int nBlock=100, int op=1, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  XPtr<BigMatrix> xpmat(pBigmat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return GenoMixer<char>(xpMat, xpmat, sirIdx, damIdx, nBlock, op, threads);
  case 2:
    return GenoMixer<short>(xpMat, xpmat, sirIdx, damIdx, nBlock, op, threads);
  case 4:
    return GenoMixer<int>(xpMat, xpmat, sirIdx, damIdx, nBlock, op, threads);
  case 8:
    return GenoMixer<double>(xpMat, xpmat, sirIdx, damIdx, nBlock, op, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}
