#include "simer.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]

template<typename T>
arma::mat calConf(XPtr<BigMatrix> pMat, int threads=0, bool verbose=true) {
  omp_setup(threads);
  
  if (verbose) { Rcout << " Computing Mendel Conflict Matrix..." << endl; }
  
  MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
  size_t i, j, k, n = pMat->nrow(), m = pMat->ncol();
  
  arma::mat numConfs(n, n, fill::zeros);

  MinimalProgressBar pb;
  Progress p(m, verbose, pb);
  
  #pragma omp parallel for schedule(dynamic) private(i, j, k)
  for (k = 0; k < m; k++) {
    for (i = 0; i < n; i++) {
      for (j = i+1; j < n; j++) {
        if (((bigm[k][i] == 0) && (bigm[k][j] == 2)) || ((bigm[k][i] == 2) && (bigm[k][k] == 0))) {
          numConfs(i, j) = numConfs(i, j) + 1;
          numConfs(i, j) = numConfs(i, j) + 1;
        }
      }
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
  StringVector kidID = rawPed[0], sirOriID = rawPed[1], damOriID = rawPed[2], sirID = rawPed[3], damID = rawPed[4], sirState = rawPed[5], damState = rawPed[6];
  size_t n = kidID.size(), m = pMat->nrow();

  StringVector fullSirID, fullDamID;
  copy(sirID.begin(), sirID.end(), back_inserter(fullSirID));
  copy(damID.begin(), damID.end(), back_inserter(fullDamID));

  StringVector candSir, candDam;
  if (candSirID.isNotNull()) {
    StringVector candSirIDUse = as<StringVector>(candSirID);
    size_t ncsu = candSirIDUse.size();
    for (size_t i = 0; i < ncsu; i++)
      fullSirID.insert(fullSirID.end(), candSirIDUse[i]);
  }
  if (candDamID.isNotNull()) {
    StringVector candDamIDUse = as<StringVector>(candDamID);
    size_t ncdu = candDamIDUse.size();
    for (size_t i = 0; i < ncdu; i++)
      fullSirID.insert(fullSirID.end(), candDamIDUse[i]);
  }
  NumericVector birdate;
  if (birthDate.isNotNull()) {
    birdate = as<NumericVector>(birthDate);
  }
  fullSirID = sort_unique(fullSirID); fullSirID.erase(0);
  fullDamID = sort_unique(fullDamID); fullDamID.erase(0);
  
  // kids should not be same as parents
  NumericVector kidOrder, sirOrder, damOrder;
  kidOrder = match(kidID, genoID); kidOrder = kidOrder - 1;
  sirOrder = match(sirID, genoID); sirOrder = sirOrder - 1;
  damOrder = match(damID, genoID); damOrder = damOrder - 1;
  LogicalVector naKid, naSir, naDam;
  naKid = is_na(kidOrder); naSir = is_na(sirOrder); naDam = is_na(damOrder);
  NumericVector sirNumConfs(n), damNumConfs(n);
  
  int exclMax = exclThres * m, assignMax = assignThres * m;

  // ******* 02 calculate conflict of pedigree in the rawPed *******
  arma::mat numConfs = calConf(pMat, threads, verbose);
  // arma::mat numConfs(pMat->ncol(), pMat->ncol(), fill::zeros);
  
  for (size_t i = 0; i < n; i++) {

    if (naKid[i]) { continue; }

    if (!naSir[i]) {
      sirNumConfs[i] = numConfs(kidOrder[i], sirOrder[i]);
      if (sirNumConfs[i] > exclMax) {
        sirState[i] = "NotFoundByGeno";
        sirID[i] = "0";
      }
    }

    if (!naDam[i]) {
      damNumConfs[i] = numConfs(kidOrder[i], damOrder[i]);
      if (damNumConfs[i] > exclMax) {
        damState[i] = "NotFoundByGeno";
        damID[i] = "0";
      }
    }

  }

  // ******* 03 seek parents of NotMatch in the rawPed *******
  StringVector candKid(n);
  LogicalVector kidFlag;
  NumericVector candKidOrder, candParOrder;
  arma::uvec candParUse;
  size_t numCand;
  arma::mat subNumConfs;

  arma::uword maxPos, rowPos, colPos;
  StringVector candPar1(1), candPar2(1);
  
  size_t i, j;

  MinimalProgressBar pb;
  Progress p(n, verbose, pb);

  if(verbose) { Rcout << " Seeking Parents..." << endl; }
  // #pragma omp parallel for schedule(dynamic) private(i, j)
  for (i = 0; i < n; i++) {

    if ((sirState[i] != "NotFoundByGeno") && (damState[i] != "NotFoundByGeno")) { continue; }
    
    candKid.fill(kidID[i]);
    kidFlag = (sirID == candKid | damID == candKid) & !naKid;
    if (birthDate.isNotNull()) {  kidFlag = kidFlag |  birdate > birdate[i]; }
    candKidOrder = kidOrder[kidFlag];
    candParOrder = wrap(arma::find(numConfs.row(kidOrder[i]) < assignMax));
    candParUse = as<arma::uvec>(setdiff(candParOrder, candKidOrder));
    numCand = candParUse.size();
    if (numCand == 0) { continue; }
    subNumConfs = numConfs.rows(candParUse);
    subNumConfs = subNumConfs.cols(candParUse);
    
    arma::uvec sortIdx = sort_index(subNumConfs);
    size_t nsi = sortIdx.n_elem;
    for (j = 0; j < nsi; j++) {
      
      maxPos = sortIdx[nsi-1-j];;
      rowPos = (maxPos + 1) % numCand;
      colPos = (maxPos + 1) / numCand;
      if (rowPos == 0) {
        rowPos = numCand;
        colPos = colPos - 1;
      }
      rowPos = rowPos - 1;
      candPar1[0] = genoID[candParUse[rowPos]];
      candPar2[0] = genoID[candParUse[colPos]];

      if ((sirState[i] == "Match") || (sirState[i] == "FoundByGeno")) {
        if (candPar1[0] != sirID[i]) {
          continue;
        } 
      }
      if ((damState[i] == "Match") || (damState[i] == "FoundByGeno")) {
        if (candPar2[0] != damID[i]) {
          continue;
        } 
      }
      
      if (find(fullSirID.begin(), fullSirID.end(), candPar1[0]) != fullSirID.end()) {
        if (find(fullDamID.begin(), fullDamID.end(), candPar2[0]) != fullDamID.end()) {
          if (sirState[i] == "NotFoundByGeno") {
            sirID[i] = candPar1[0];
            sirState[i] = "FoundByGeno";
            sirNumConfs[i] = numConfs(kidOrder[i], candParUse[rowPos]);
          }
          if (damState[i] == "NotFoundByGeno") {
            damID[i] = candPar2[0];
            damState[i] = "FoundByGeno";
            damNumConfs[i] = numConfs(kidOrder[i], candParUse[colPos]);
          }
          if (((sirState[i] == "Match") || (sirState[i] == "FoundByGeno")) && ((damState[i] == "Match") || (damState[i] == "FoundByGeno"))) {
            break;
          }
        }
      }
      
    }

    if ( ! Progress::check_abort() ) { p.increment(); }
  }
  
  DataFrame parConflict = DataFrame::create(
    Named("kid")        = kidID,
    _["sirOrigin"]      = sirOriID,
    _["damOrigin"]      = damOriID,
    _["sirFound"]       = sirID,
    _["damFound"]       = damID,
    _["sirState"]       = sirState,
    _["damState"]       = damState,
    _["sirNumConfs"]    = sirNumConfs,
    _["damNumConfs"]    = damNumConfs,
    _["sirRatioConfs"]  = sirNumConfs * 100 / m,
    _["damRatioConfs"]  = damNumConfs * 100 / m
  );
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

