#include <RcppArmadillo.h>
#include "simer_omp.h"
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <R_ext/Print.h>
#include <progress.hpp>
#include "progress_bar.hpp"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace std;
using namespace Rcpp;
using namespace arma;

class MinimalProgressBar: public ProgressBar{
  public:
  MinimalProgressBar()  {
    _finalized = false;
  }
  ~MinimalProgressBar() {}
  void display() {}
  void update(float progress) {
    if (_finalized) return;
    REprintf("\r");
    REprintf(" Calculating in process...(finished %.2f%)", progress * 100);
  }
  void end_display() {
    if (_finalized) return;
    REprintf("\r");
    
    REprintf(" Calculating in process...(finished 100.00%)");
    REprintf("\n");
    _finalized = true;
  }
  private:
  bool _finalized;
};

template<typename T>
arma::mat calConf(XPtr<BigMatrix> pMat, int threads=0, bool verbose=true) {
  omp_setup(threads);
  
  if(verbose) { Rcout << " Computing Mendel Conflict Matrix..." << endl; }
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  size_t m = pMat->nrow();
  size_t n = pMat->ncol();
  size_t i, j, k;
  
  arma::mat numConfs(n, n, fill::zeros);
  arma::vec coli(m), colj(m);
  
  MinimalProgressBar pb;
  Progress p(n, verbose, pb);
  
  #pragma omp parallel for schedule(dynamic) firstprivate(coli, colj) private(i, j, k)
  for (i = 0; i < n; i++) {
    for(k = 0; k < m; k++){
      coli[k] = bigm[i][k];
    }
    coli.elem(find(coli == 1)).fill(3);
    for (j = i+1; j < n; j++) {
      for (k = 0; k < m; k++) {
        colj[k] = bigm[j][k];
      }
      arma::uvec confIdx = arma::find(coli == 2 - colj);
      numConfs(i, j) = confIdx.size();
      numConfs(j, i) = confIdx.size();
    }
    if ( ! Progress::check_abort() ) { p.increment(); }
  }
  
  return numConfs;
}

arma::mat calConf(XPtr<BigMatrix> pMat, int threads=0, bool verbose=true) {
  XPtr<BigMatrix> xpMat(pMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return calConf<char>(xpMat, threads, verbose);
  case 2:
    return calConf<short>(xpMat, threads, verbose);
  case 4:
    return calConf<int>(xpMat, threads, verbose);
  case 8:
    return calConf<double>(xpMat, threads, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template <typename T>
DataFrame PedigreeCorrector(XPtr<BigMatrix> pMat, StringVector genoID, DataFrame rawPed, Nullable<StringVector> candSirID=R_NilValue, Nullable<StringVector> candDamID=R_NilValue, double exclThres=0.005, double assignThres=0.01, Nullable<NumericVector> birthDate=R_NilValue, int threads=0, bool verbose=true) {
  omp_setup(threads);
  
  // ******* 01 prepare data for checking rawPed *******
  StringVector kidID = rawPed[0], sirID = rawPed[1], damID = rawPed[2];
  size_t n = kidID.size(), m = pMat->nrow();
  NumericVector kidOrder, sirOrder, damOrder;
  kidOrder = match(kidID, genoID); kidOrder = kidOrder - 1;
  sirOrder = match(sirID, genoID); sirOrder = sirOrder - 1; 
  damOrder = match(damID, genoID); damOrder = damOrder - 1;
  LogicalVector naKid, naSir, naDam;
  naKid = is_na(kidOrder); naSir = is_na(sirOrder); naDam = is_na(damOrder);
  StringVector sirState(n), damState(n);
  NumericVector sirNumConfs(n), damNumConfs(n);
  
  int exclMax = exclThres * m, assignMax = assignThres * m;
  
  // ******* 02 check rawPed *******
  StringVector zero(n, "0");
  sirState[naSir] = "NotFound"; damState[naDam] = "NotFound";
  sirState[naKid] = "NoGeno"; damState[naKid] = "NoGeno";
  sirState[(sirID == zero) & (damID == zero)] = "NoGeno"; 
  damState[(sirID == zero) & (damID == zero)] = "NoGeno"; 

  // kids should not be same as parents
  sirState[kidID == sirID] = "NotFound";
  damState[kidID == damID] = "NotFound";

  // calculate conflict of pedigree in the rawPed
  arma::mat numConfs = calConf(pMat, threads, verbose);
  // arma::mat numConfs(pMat->ncol(), pMat->ncol(), fill::zeros);
  
  for (size_t i = 0; i < n; i++) {

    if (naKid[i]) { continue; }

    if (!naSir[i]) {
      sirNumConfs[i] = numConfs(kidOrder[i], sirOrder[i]);
      if (sirNumConfs[i] <= exclMax) {
        sirState[i] = "Match";
      } else {
        sirState[i] = "NotFound";
      }
    } // if (!isnan(sirOrder[i])) {

    if (!naDam[i]) {
      damNumConfs[i] = numConfs(kidOrder[i], damOrder[i]);
      if (damNumConfs[i] <= exclMax) {
        damState[i] = "Match";
      } else {
        damState[i] = "NotFound";
      }
    } // if (!isnan(damOrder[i])) {

  } // for (int i = 0; i < n; i++) {

  // ******* 03 prepare data for seeking parents *******
  StringVector fullSirID, fullDamID;
  copy(sirID.begin(), sirID.end(), back_inserter(fullSirID));
  copy(damID.begin(), damID.end(), back_inserter(fullDamID));

  StringVector candSir, candDam;
  if (candSirID.isNotNull()) {
    StringVector candSirIDUse = as<StringVector>(candSirID);
    for (size_t i = 0; i < candSirIDUse.size(); i++)
      fullSirID.insert(fullSirID.end(), candSirIDUse[i]);
  }
  if (candDamID.isNotNull()) {
    StringVector candDamIDUse = as<StringVector>(candDamID);
    for (size_t i = 0; i < candDamIDUse.size(); i++)
      fullSirID.insert(fullSirID.end(), candDamIDUse[i]);
  }
  NumericVector birdate;
  if (birthDate.isNotNull()) {
    birdate = as<NumericVector>(birthDate);
  }

  fullSirID = sort_unique(fullSirID);
  fullDamID = sort_unique(fullDamID);
  
  StringVector candKid(n);
  LogicalVector kidFlag;
  NumericVector candKidOrder, candParOrder;
  arma::uvec candParUse;
  size_t numCand;
  arma::mat subNumConfs;

  arma::uvec findPos;
  arma::uword maxPos, rowPos, colPos;
  string candPar1, candPar2;
  
  size_t i, j;

  MinimalProgressBar pb;
  Progress p(n, verbose, pb);

  // ******* 04 seek parents of NotMatch in the rawPed *******
  if(verbose) { Rcout << " Seeking Parents..." << endl; }
  // #pragma omp parallel for schedule(dynamic) private(i, j)
  for (i = 0; i < n; i++) {

    if ((sirState[i] != "NotFound") && (damState[i] != "NotFound")) { continue; }
    
    candKid.fill(kidID[i]);
    kidFlag = (sirID == candKid | damID == candKid) & !naKid;
    if (birthDate.isNotNull()) {  kidFlag = kidFlag |  birdate > birdate[i]; }
    candKidOrder = kidOrder[kidFlag];
    candParOrder = wrap(arma::find(numConfs.row(kidOrder[i]) < assignMax));
    candParUse = as<arma::uvec>(setdiff(candParOrder, candKidOrder));
    numCand = candParUse.size();
    if (numCand == 0) { continue;  }
    subNumConfs = numConfs.rows(candParUse);
    subNumConfs = subNumConfs.cols(candParUse);

    arma::uvec sortIdx = sort_index(subNumConfs);
    for (j = sortIdx.max(); j > 0; j--) {
      
      findPos = arma::find(sortIdx == j);
      maxPos = findPos[0];
      rowPos = (maxPos + 1) % numCand;
      colPos = (maxPos + 1) / numCand;
      if (rowPos == 0) {
        rowPos = numCand;
        colPos = colPos - 1;
      }
      rowPos = rowPos - 1;
      candPar1 = as<string>(genoID[candParUse[rowPos]]);
      candPar2 = as<string>(genoID[candParUse[colPos]]);

      if (find(fullSirID.begin(), fullSirID.end(), candPar1) != fullSirID.end()) {
        if (find(fullDamID.begin(), fullDamID.end(), candPar2) != fullDamID.end()) {
          if (candPar1.compare(sirID[i])) {
            sirID[i] = candPar1;
            sirState[i] = "Found";
            sirNumConfs[i] = numConfs(kidOrder[i], candParUse[rowPos]);
          }
          if (candPar2.compare(damID[i])) {
            damID[i] = candPar2;
            damState[i] = "Found";
            damNumConfs[i] = numConfs(kidOrder[i], candParUse[colPos]);
          }
          break;
        }

      } else if (find(fullSirID.begin(), fullSirID.end(), candPar2) != fullSirID.end()) {
        if (find(fullDamID.begin(), fullDamID.end(), candPar1) != fullDamID.end()) {
          if (candPar2.compare(sirID[i])) {
            sirID[i] = candPar2;
            sirState[i] = "Found";
            sirNumConfs[i] = numConfs(kidOrder[i], candParUse[colPos]);
          }
          if (candPar1.compare(damID[i])) {
            damID[i] = candPar1;
            damState[i] = "Found";
            damNumConfs[i] = numConfs(kidOrder[i], candParUse[rowPos]);
          }
          break;
        }
      }
      
    }

    if ( ! Progress::check_abort() ) { p.increment(); }
  }
  
  DataFrame parConflict = DataFrame::create(Named("kid") = kidID,
                                                _["sir"] = sirID, 
                                                _["dam"] = damID, 
                                                _["sirState"] = sirState, 
                                                _["damState"] = damState, 
                                                _["sirNumConfs"] = sirNumConfs, 
                                                _["damNumConfs"] = damNumConfs, 
                                                _["sirRatioConfs"] = sirNumConfs * 100 / m, 
                                                _["damRatioConfs"] = damNumConfs * 100 / m);
  return parConflict;
}

// [[Rcpp::export]]
DataFrame PedigreeCorrector(const SEXP pBigMat, StringVector rawGenoID, DataFrame rawPed, Nullable<StringVector> candSirID=R_NilValue, Nullable<StringVector> candDamID=R_NilValue, double exclThres=0.005, double assignThres=0.02, Nullable<NumericVector> birthDate=R_NilValue, int threads=0, bool verbose=true){
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return PedigreeCorrector<char>(xpMat, rawGenoID, rawPed, candSirID, candDamID, exclThres, assignThres, birthDate, threads, verbose);
  case 2:
    return PedigreeCorrector<short>(xpMat, rawGenoID, rawPed, candSirID, candDamID, exclThres, assignThres, birthDate, threads, verbose);
  case 4:
    return PedigreeCorrector<int>(xpMat, rawGenoID, rawPed, candSirID, candDamID, exclThres, assignThres, birthDate, threads, verbose);
  case 8:
    return PedigreeCorrector<double>(xpMat, rawGenoID, rawPed, candSirID, candDamID, exclThres, assignThres, birthDate, threads, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

