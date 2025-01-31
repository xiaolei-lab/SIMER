#include "simer.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

template<typename T>
NumericVector FilterMind(XPtr<BigMatrix> pMat, double NA_C, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  
  size_t i, j, n = pMat->nrow(), m = pMat->ncol();
  NumericVector indNumNA(n, 0);
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      if (bigm[j][i] == NA_C) { indNumNA[i] += 1;  }
    }
  }
  
  return indNumNA;
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
NumericVector FilterGeno(XPtr<BigMatrix> pMat, double NA_C, IntegerVector &indIdx, IntegerVector &mrkIdx, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  size_t i, j, n = indIdx.size(), m = mrkIdx.size();
  NumericVector mrkNumNA(m, 0);
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      if (bigm[mrkIdx[j]][indIdx[i]] == NA_C) { mrkNumNA[j] += 1;  }
    }
  }
  
  return mrkNumNA;
}

NumericVector FilterGeno(const SEXP pBigMat, double NA_C, IntegerVector &indIdx, IntegerVector &mrkIdx, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return FilterGeno<char>(xpMat, NA_CHAR, indIdx, mrkIdx, threads);
  case 2:
    return FilterGeno<short>(xpMat, NA_SHORT, indIdx, mrkIdx, threads);
  case 4:
    return FilterGeno<int>(xpMat, NA_INTEGER, indIdx, mrkIdx, threads);
  case 8:
    return FilterGeno<double>(xpMat, NA_REAL, indIdx, mrkIdx, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
arma::mat CalGenoFreq(XPtr<BigMatrix> pMat, double NA_C, IntegerVector &indIdx, IntegerVector &mrkIdx, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  
  size_t i, j, n = indIdx.size(), m = mrkIdx.size();
  arma::mat genoFreq(m, 3, fill::zeros);

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      if (bigm[mrkIdx[j]][indIdx[i]] == 0) {
        genoFreq(j, 0) = genoFreq(j, 0) + 1; 
      } else if (bigm[mrkIdx[j]][indIdx[i]] == 1) {
        genoFreq(j, 1) = genoFreq(j, 1) + 1;
      } else if (bigm[mrkIdx[j]][indIdx[i]] == 2) {
        genoFreq(j, 2) = genoFreq(j, 2) + 1;
      }
    }
  }
  
  return genoFreq;
}

arma::mat CalGenoFreq(const SEXP pBigMat, IntegerVector &indIdx, IntegerVector &mrkIdx, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return CalGenoFreq<char>(xpMat, NA_CHAR, indIdx, mrkIdx, threads);
  case 2:
    return CalGenoFreq<short>(xpMat, NA_SHORT, indIdx, mrkIdx, threads);
  case 4:
    return CalGenoFreq<int>(xpMat, NA_INTEGER, indIdx, mrkIdx, threads);
  case 8:
    return CalGenoFreq<double>(xpMat, NA_REAL, indIdx, mrkIdx, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

NumericVector FilterMAF(arma::mat &genoFreq, int threads=0) {
  omp_setup(threads);
  
  IntegerVector freq0 = wrap(genoFreq.col(0));
  IntegerVector freq1 = wrap(genoFreq.col(1));
  IntegerVector freq2 = wrap(genoFreq.col(2));

  size_t i, m = genoFreq.n_rows;
  NumericVector MAF(m); MAF.fill(0);
  
  #pragma omp parallel for schedule(dynamic) private(i)
  for (i = 0; i < m; i++) {
    MAF[i] = (freq0[i] + freq1[i] * 0.5) / 
      (freq0[i] + freq1[i] + freq2[i]);
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

NumericVector FilterHWE(arma::mat &genoFreq, int threads=0) {
  omp_setup(threads);
  
  size_t i, m = genoFreq.n_rows;
  IntegerVector freq0 = wrap(genoFreq.col(0));
  IntegerVector freq1 = wrap(genoFreq.col(1));
  IntegerVector freq2 = wrap(genoFreq.col(2));
  NumericVector PVAL(m); PVAL.fill(0);
  
  #pragma omp parallel for schedule(dynamic) private(i)
  for (i = 0; i < m; i++) {
    PVAL[i] = SNPHWE(freq1[i], freq0[i], freq2[i]);
  }
  
  return PVAL;
}

// [[Rcpp::export]]
List GenoFilter(const SEXP pBigMat, double NA_C, Nullable<IntegerVector> keepIndsNull=R_NilValue, Nullable<double> filterGeno=R_NilValue, Nullable<double> filterHWE=R_NilValue, Nullable<double> filterMind=R_NilValue, Nullable<double> filterMAF=R_NilValue, int threads=0, bool verbose=true) {
  omp_setup(threads);

  XPtr<BigMatrix> pMat(pBigMat);

  double n = pMat->nrow(), m = pMat->ncol();
  IntegerVector  keepMrks = seq(0, m - 1);
  IntegerVector keepInds;
  if (keepIndsNull.isNull()) {
    keepInds = seq(0, n - 1);
  } else {
    keepInds = as<IntegerVector>(keepIndsNull);
    keepInds = keepInds - 1;
    n = keepInds.size();
  }
  
  double fgeno = 0, fhwe = 0, fmaf = 0, fmind = 0;
  if (filterGeno.isNotNull()) { fgeno = as<double>(filterGeno); }
  if (filterHWE.isNotNull() ) { fhwe  = as<double>(filterHWE ); }
  if (filterMAF.isNotNull() ) { fmaf  = as<double>(filterMAF); }
  if (filterMind.isNotNull()) { fmind = as<double>(filterMind); }
  
  if (verbose) {
    Rcout << " Options in effect:" << endl;
    if (keepIndsNull.isNotNull()  ) { Rcout << "   --keep-ind filePed "       << endl; }
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
    NumericVector indNumNA = FilterMind(pBigMat, NA_C, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepInds = keepInds[indNumNA/m < fmind];
    if (verbose) {
      Rcout << " " << (n - keepInds.size())  << " samples removed due to missing genotype data (--mind)." << endl;
      n = keepInds.size();
      Rcout << " " << n << " samples remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  if (filterGeno.isNotNull()) {
    if (verbose) { Rcout << " Calculating variant missingness rates..."; }
    NumericVector mrkNumNA = FilterGeno(pBigMat, NA_C, keepInds, keepMrks, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepMrks = keepMrks[mrkNumNA/n < fgeno];
    if (verbose) {
      Rcout << " " << (m - keepMrks.size()) << " variants removed due to missing genotype data (--geno)." << endl;
      m = keepMrks.size();
      Rcout << " " << m << " variants remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  arma::mat genoFreq;
  if (filterMAF.isNotNull() || filterHWE.isNotNull()) {
    if (verbose) { Rcout << " Calculating Genotype FrequenindIdxes..."; }
    genoFreq = CalGenoFreq(pBigMat, keepInds, keepMrks, threads);
    if (verbose) {  Rcout << " done." << endl << endl; }
  }
  
  if (filterHWE.isNotNull()) {
    if (verbose) { Rcout << " Performing Hardy-Weinberg test..."; }
    NumericVector PVAL = FilterHWE(genoFreq, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepMrks = keepMrks[PVAL > fhwe];
    arma::vec armaPVAL = as<arma::vec>(PVAL);
    genoFreq = genoFreq.rows(arma::find(armaPVAL > fhwe));
    if (verbose) {
      Rcout << " " << (m - keepMrks.size()) << " variants removed due to exceeding HWE-P-Value (--hwe)." << endl;
      m = keepMrks.size();
      Rcout << " " << m << " variants remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  if (filterMAF.isNotNull()) {
    if (verbose) { Rcout << " Calculating Minor Allele FrequenindIdxes..."; }
    NumericVector MAF = FilterMAF(genoFreq, threads);
    if (verbose) {  Rcout << " done." << endl; }
    keepMrks = keepMrks[MAF >= fmaf];
    if (verbose) {
      Rcout << " " << (m - keepMrks.size()) << " variants removed due to exceeding MAF (--maf)." << endl;
      m = keepMrks.size();
      Rcout << " " << m << " variants remaining after main filters." << endl;
      Rcout << endl;
    }
  }
  
  keepMrks = keepMrks + 1;
  keepInds = keepInds + 1;
  List genoInfo = List::create(Named("keepMrks") = keepMrks,
                                   _["keepInds"] = keepInds);
  return genoInfo;
}

template<typename T>
void Mat2BigMat(XPtr<BigMatrix> pMat, IntegerMatrix &mat, Nullable<IntegerVector> indIdxNull=R_NilValue, int op=1, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  
  IntegerVector indIdx;
  if (indIdxNull.isNull()) {
    indIdx = seq(0, mat.nrow() - 1);
  } else {
    indIdx = as<IntegerVector>(indIdxNull);
    indIdx = indIdx - 1;
  }
  
  size_t i, j, n = indIdx.length(), m = mat.ncol(), n0 = pMat->nrow(), m0 = pMat->ncol();
  op = op - 1;
  if (m != m0) {
    Rcpp::stop("'bigmat' and 'mat' should have the same marker number!");
  }
  if (op + n > n0) {
    Rcpp::stop("'mat' cannot be intert to bigmat completely!");
  }
  if (max(indIdx) + 1 > mat.nrow()) {
    Rcpp::stop("'indIdx' is out of bound!");
  }
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      bigmat[j][op + i] = mat(indIdx[i], j);
    }
  }

}

// [[Rcpp::export]]
void Mat2BigMat(const SEXP pBigMat, IntegerMatrix &mat, Nullable<IntegerVector> indIdxNull=R_NilValue, int op=1, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return Mat2BigMat<char>(xpMat, mat, indIdxNull, op, threads);
  case 2:
    return Mat2BigMat<short>(xpMat, mat, indIdxNull, op, threads);
  case 4:
    return Mat2BigMat<int>(xpMat, mat, indIdxNull, op, threads);
  case 8:
    return Mat2BigMat<double>(xpMat, mat, indIdxNull, op, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void BigMat2BigMat(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pmat, Nullable<IntegerVector> indIdxNull=R_NilValue, int op=1, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pmat);
  
  IntegerVector indIdx;
  size_t n1 = pmat->nrow();
  if (indIdxNull.isNull()) {
    indIdx = seq(0, n1 - 1);
  } else {
    indIdx = as<IntegerVector>(indIdxNull);
    indIdx = indIdx - 1;
  }
  
  size_t i, j, n = indIdx.length(), m = pmat->ncol(), n0 = pMat->nrow(), m0 = pMat->ncol();
  op = op - 1;
  if (m != m0) {
    Rcpp::stop("'bigmat' and 'pmat' should have the same marker number!");
  }
  if (op + n > n0) {
    Rcpp::stop("'pmat' cannot be intert to bigmat completely!");
  }
  size_t maxIndIdx = max(indIdx);
  if (maxIndIdx + 1 > n1) {
    Rcpp::stop("'indIdx' is out of bound!");
  }
  
  IntegerMatrix mat(n1, m);
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n1; i++) {
      mat(i, j) = bigm[j][i];
    }
  }
  
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      bigmat[j][op + i] = mat(indIdx[i], j);
    }
  }
  
}

// [[Rcpp::export]]
void BigMat2BigMat(const SEXP pBigMat, const SEXP pBigmat, Nullable<IntegerVector> indIdxNull=R_NilValue, int op=1, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  XPtr<BigMatrix> xpmat(pBigmat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return BigMat2BigMat<char>(xpMat, xpmat, indIdxNull, op, threads);
  case 2:
    return BigMat2BigMat<short>(xpMat, xpmat, indIdxNull, op, threads);
  case 4:
    return BigMat2BigMat<int>(xpMat, xpmat, indIdxNull, op, threads);
  case 8:
    return BigMat2BigMat<double>(xpMat, xpmat, indIdxNull, op, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void geno_cvt1_mat(XPtr<BigMatrix> pMat, IntegerMatrix &mat, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);

  size_t i, j, n = pMat->nrow(), m = pMat->ncol();

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      bigmat[j][i] = mat(2*i, j) + mat(2*i+1, j);
    }
  }
  
}

// [[Rcpp::export]]
void geno_cvt1_mat(const SEXP pBigMat, IntegerMatrix &mat, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return geno_cvt1_mat<char>(xpMat, mat, threads);
  case 2:
    return geno_cvt1_mat<short>(xpMat, mat, threads);
  case 4:
    return geno_cvt1_mat<int>(xpMat, mat, threads);
  case 8:
    return geno_cvt1_mat<double>(xpMat, mat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void geno_cvt1_bigmat(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pmat, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pmat);

  size_t i, j, n = pMat->nrow(), m = pMat->ncol();

  IntegerMatrix mat(n, m);
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      mat(i, j) = bigm[j][2*i] + bigm[j][2*i+1];
    }
  }

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      bigmat[j][i] = mat(i, j);
    }
  }
  
}

// [[Rcpp::export]]
void geno_cvt1_bigmat(const SEXP pBigMat, const SEXP pBigmat, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  XPtr<BigMatrix> xpmat(pBigmat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return geno_cvt1_bigmat<char>(xpMat, xpmat, threads);
  case 2:
    return geno_cvt1_bigmat<short>(xpMat, xpmat, threads);
  case 4:
    return geno_cvt1_bigmat<int>(xpMat, xpmat, threads);
  case 8:
    return geno_cvt1_bigmat<double>(xpMat, xpmat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void geno_cvt2_mat(XPtr<BigMatrix> pMat, IntegerMatrix &mat, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);

  size_t i, j, n = mat.nrow(), m = mat.ncol();

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      if (mat(i, j) == 0) {
        bigmat[j][2*i] = 0;
        bigmat[j][2*i+1] = 0;
      } else if (mat(i, j) == 1) {
        bigmat[j][2*i] = 0;
        bigmat[j][2*i+1] = 1;
      } else if (mat(i, j) == 2) {
        bigmat[j][2*i] = 1;
        bigmat[j][2*i+1] = 1;
      } else {
        Rcpp::stop("Elements in genotype data should be 0, 1 or 2!");
      }
    }
  }
  
}

// [[Rcpp::export]]
void geno_cvt2_mat(const SEXP pBigMat, IntegerMatrix &mat, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return geno_cvt2_mat<char>(xpMat, mat, threads);
  case 2:
    return geno_cvt2_mat<short>(xpMat, mat, threads);
  case 4:
    return geno_cvt2_mat<int>(xpMat, mat, threads);
  case 8:
    return geno_cvt2_mat<double>(xpMat, mat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void geno_cvt2_bigmat(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pmat, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pmat);

  size_t i, j, n = pmat->nrow(), m = pmat->ncol();

  IntegerMatrix mat(2 * n, m);
  #pragma omp parallel for schedule(dynamic) private(i, j)
    for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      if (bigm[j][i] == 0) {
        mat(2*i, j) = 0;
        mat(2*i+1, j) = 0;
      } else if (bigm[j][i] == 1) {
        mat(2*i, j) = 0;
        mat(2*i+1, j) = 1;
      } else if (bigm[j][i] == 2) {
        mat(2*i, j) = 1;
        mat(2*i+1, j) = 1;
      } else {
        Rcpp::stop("Elements in genotype data should be 0, 1 or 2!");
      }
    }
  }

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < 2 * n; i++) {
      bigmat[j][i] = mat(i, j);
    }
  }
  
}

// [[Rcpp::export]]
void geno_cvt2_bigmat(const SEXP pBigMat, const SEXP pBigmat, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  XPtr<BigMatrix> xpmat(pBigmat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return geno_cvt2_bigmat<char>(xpMat, xpmat, threads);
  case 2:
    return geno_cvt2_bigmat<short>(xpMat, xpmat, threads);
  case 4:
    return geno_cvt2_bigmat<int>(xpMat, xpmat, threads);
  case 8:
    return geno_cvt2_bigmat<double>(xpMat, xpmat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void bigt_mat(XPtr<BigMatrix> pMat, IntegerMatrix &mat, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);

  size_t i, j, n = pMat->nrow(), m = pMat->ncol();

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      bigmat[j][i] = mat(j, i);
    }
  }
  
}

// [[Rcpp::export]]
void bigt_mat(const SEXP pBigMat, IntegerMatrix &mat, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return bigt_mat<char>(xpMat, mat, threads);
  case 2:
    return bigt_mat<short>(xpMat, mat, threads);
  case 4:
    return bigt_mat<int>(xpMat, mat, threads);
  case 8:
    return bigt_mat<double>(xpMat, mat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template<typename T>
void bigt_bigmat(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pmat, int threads=0) {
  omp_setup(threads);
  
  MatrixAccessor<T> bigmat = MatrixAccessor<T>(*pMat);
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pmat);

  size_t i, j, n = pMat->nrow(), m = pMat->ncol();

  IntegerMatrix mat(2 * n, m);
  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      mat(i, j) = bigm[i][j];
    }
  }

  #pragma omp parallel for schedule(dynamic) private(i, j)
  for (j = 0; j < m; j++) {
    for (i = 0; i < n; i++) {
      bigmat[j][i] = mat(i, j);
    }
  }
  
}

// [[Rcpp::export]]
void bigt_bigmat(const SEXP pBigMat, const SEXP pBigmat, int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  XPtr<BigMatrix> xpmat(pBigmat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return bigt_bigmat<char>(xpMat, xpmat, threads);
  case 2:
    return bigt_bigmat<short>(xpMat, xpmat, threads);
  case 4:
    return bigt_bigmat<int>(xpMat, xpmat, threads);
  case 8:
    return bigt_bigmat<double>(xpMat, xpmat, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}
