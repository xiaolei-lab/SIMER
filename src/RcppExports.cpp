// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// write_bfile
void write_bfile(SEXP pBigMat, std::string bed_file, bool mrkbycol, int threads, bool verbose);
RcppExport SEXP _simer_write_bfile(SEXP pBigMatSEXP, SEXP bed_fileSEXP, SEXP mrkbycolSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< std::string >::type bed_file(bed_fileSEXP);
    Rcpp::traits::input_parameter< bool >::type mrkbycol(mrkbycolSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    write_bfile(pBigMat, bed_file, mrkbycol, threads, verbose);
    return R_NilValue;
END_RCPP
}
// read_bfile
void read_bfile(std::string bed_file, SEXP pBigMat, long maxLine, int threads, bool verbose);
RcppExport SEXP _simer_read_bfile(SEXP bed_fileSEXP, SEXP pBigMatSEXP, SEXP maxLineSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bed_file(bed_fileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< long >::type maxLine(maxLineSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    read_bfile(bed_file, pBigMat, maxLine, threads, verbose);
    return R_NilValue;
END_RCPP
}
// GenoFilter
List GenoFilter(const SEXP pBigMat, double NA_C, Nullable<IntegerVector> keepIndsNull, Nullable<double> filterGeno, Nullable<double> filterHWE, Nullable<double> filterMind, Nullable<double> filterMAF, int threads, bool verbose);
RcppExport SEXP _simer_GenoFilter(SEXP pBigMatSEXP, SEXP NA_CSEXP, SEXP keepIndsNullSEXP, SEXP filterGenoSEXP, SEXP filterHWESEXP, SEXP filterMindSEXP, SEXP filterMAFSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< double >::type NA_C(NA_CSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type keepIndsNull(keepIndsNullSEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type filterGeno(filterGenoSEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type filterHWE(filterHWESEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type filterMind(filterMindSEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type filterMAF(filterMAFSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(GenoFilter(pBigMat, NA_C, keepIndsNull, filterGeno, filterHWE, filterMind, filterMAF, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// Mat2BigMat
void Mat2BigMat(const SEXP pBigMat, IntegerMatrix& mat, Nullable<IntegerVector> indIdxNull, int op, int threads);
RcppExport SEXP _simer_Mat2BigMat(SEXP pBigMatSEXP, SEXP matSEXP, SEXP indIdxNullSEXP, SEXP opSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type indIdxNull(indIdxNullSEXP);
    Rcpp::traits::input_parameter< int >::type op(opSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Mat2BigMat(pBigMat, mat, indIdxNull, op, threads);
    return R_NilValue;
END_RCPP
}
// BigMat2BigMat
void BigMat2BigMat(const SEXP pBigMat, const SEXP pBigmat, Nullable<IntegerVector> indIdxNull, int op, int threads);
RcppExport SEXP _simer_BigMat2BigMat(SEXP pBigMatSEXP, SEXP pBigmatSEXP, SEXP indIdxNullSEXP, SEXP opSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type pBigmat(pBigmatSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type indIdxNull(indIdxNullSEXP);
    Rcpp::traits::input_parameter< int >::type op(opSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    BigMat2BigMat(pBigMat, pBigmat, indIdxNull, op, threads);
    return R_NilValue;
END_RCPP
}
// geno_cvt1_mat
void geno_cvt1_mat(const SEXP pBigMat, IntegerMatrix& mat, int threads);
RcppExport SEXP _simer_geno_cvt1_mat(SEXP pBigMatSEXP, SEXP matSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    geno_cvt1_mat(pBigMat, mat, threads);
    return R_NilValue;
END_RCPP
}
// geno_cvt1_bigmat
void geno_cvt1_bigmat(const SEXP pBigMat, const SEXP pBigmat, int threads);
RcppExport SEXP _simer_geno_cvt1_bigmat(SEXP pBigMatSEXP, SEXP pBigmatSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type pBigmat(pBigmatSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    geno_cvt1_bigmat(pBigMat, pBigmat, threads);
    return R_NilValue;
END_RCPP
}
// geno_cvt2_mat
void geno_cvt2_mat(const SEXP pBigMat, IntegerMatrix& mat, int threads);
RcppExport SEXP _simer_geno_cvt2_mat(SEXP pBigMatSEXP, SEXP matSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    geno_cvt2_mat(pBigMat, mat, threads);
    return R_NilValue;
END_RCPP
}
// geno_cvt2_bigmat
void geno_cvt2_bigmat(const SEXP pBigMat, const SEXP pBigmat, int threads);
RcppExport SEXP _simer_geno_cvt2_bigmat(SEXP pBigMatSEXP, SEXP pBigmatSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type pBigmat(pBigmatSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    geno_cvt2_bigmat(pBigMat, pBigmat, threads);
    return R_NilValue;
END_RCPP
}
// bigt_mat
void bigt_mat(const SEXP pBigMat, IntegerMatrix& mat, int threads);
RcppExport SEXP _simer_bigt_mat(SEXP pBigMatSEXP, SEXP matSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    bigt_mat(pBigMat, mat, threads);
    return R_NilValue;
END_RCPP
}
// bigt_bigmat
void bigt_bigmat(const SEXP pBigMat, const SEXP pBigmat, int threads);
RcppExport SEXP _simer_bigt_bigmat(SEXP pBigMatSEXP, SEXP pBigmatSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type pBigmat(pBigmatSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    bigt_bigmat(pBigMat, pBigmat, threads);
    return R_NilValue;
END_RCPP
}
// impute_marker
void impute_marker(SEXP pBigMat, bool mrkbycol, int threads, bool verbose);
RcppExport SEXP _simer_impute_marker(SEXP pBigMatSEXP, SEXP mrkbycolSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< bool >::type mrkbycol(mrkbycolSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    impute_marker(pBigMat, mrkbycol, threads, verbose);
    return R_NilValue;
END_RCPP
}
// hasNA
bool hasNA(SEXP pBigMat, bool mrkbycol, const Nullable<arma::uvec> geno_ind, const Nullable<arma::uvec> marker_ind, const int threads);
RcppExport SEXP _simer_hasNA(SEXP pBigMatSEXP, SEXP mrkbycolSEXP, SEXP geno_indSEXP, SEXP marker_indSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< bool >::type mrkbycol(mrkbycolSEXP);
    Rcpp::traits::input_parameter< const Nullable<arma::uvec> >::type geno_ind(geno_indSEXP);
    Rcpp::traits::input_parameter< const Nullable<arma::uvec> >::type marker_ind(marker_indSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(hasNA(pBigMat, mrkbycol, geno_ind, marker_ind, threads));
    return rcpp_result_gen;
END_RCPP
}
// hasNABed
bool hasNABed(std::string bed_file, int ind, long maxLine, int threads, bool verbose);
RcppExport SEXP _simer_hasNABed(SEXP bed_fileSEXP, SEXP indSEXP, SEXP maxLineSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bed_file(bed_fileSEXP);
    Rcpp::traits::input_parameter< int >::type ind(indSEXP);
    Rcpp::traits::input_parameter< long >::type maxLine(maxLineSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(hasNABed(bed_file, ind, maxLine, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// PedigreeCorrector
DataFrame PedigreeCorrector(const SEXP pBigMat, StringVector rawGenoID, DataFrame rawPed, Nullable<StringVector> candSirID, Nullable<StringVector> candDamID, double exclThres, double assignThres, Nullable<NumericVector> birthDate, int threads, bool verbose);
RcppExport SEXP _simer_PedigreeCorrector(SEXP pBigMatSEXP, SEXP rawGenoIDSEXP, SEXP rawPedSEXP, SEXP candSirIDSEXP, SEXP candDamIDSEXP, SEXP exclThresSEXP, SEXP assignThresSEXP, SEXP birthDateSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< StringVector >::type rawGenoID(rawGenoIDSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type rawPed(rawPedSEXP);
    Rcpp::traits::input_parameter< Nullable<StringVector> >::type candSirID(candSirIDSEXP);
    Rcpp::traits::input_parameter< Nullable<StringVector> >::type candDamID(candDamIDSEXP);
    Rcpp::traits::input_parameter< double >::type exclThres(exclThresSEXP);
    Rcpp::traits::input_parameter< double >::type assignThres(assignThresSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type birthDate(birthDateSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(PedigreeCorrector(pBigMat, rawGenoID, rawPed, candSirID, candDamID, exclThres, assignThres, birthDate, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_simer_write_bfile", (DL_FUNC) &_simer_write_bfile, 5},
    {"_simer_read_bfile", (DL_FUNC) &_simer_read_bfile, 5},
    {"_simer_GenoFilter", (DL_FUNC) &_simer_GenoFilter, 9},
    {"_simer_Mat2BigMat", (DL_FUNC) &_simer_Mat2BigMat, 5},
    {"_simer_BigMat2BigMat", (DL_FUNC) &_simer_BigMat2BigMat, 5},
    {"_simer_geno_cvt1_mat", (DL_FUNC) &_simer_geno_cvt1_mat, 3},
    {"_simer_geno_cvt1_bigmat", (DL_FUNC) &_simer_geno_cvt1_bigmat, 3},
    {"_simer_geno_cvt2_mat", (DL_FUNC) &_simer_geno_cvt2_mat, 3},
    {"_simer_geno_cvt2_bigmat", (DL_FUNC) &_simer_geno_cvt2_bigmat, 3},
    {"_simer_bigt_mat", (DL_FUNC) &_simer_bigt_mat, 3},
    {"_simer_bigt_bigmat", (DL_FUNC) &_simer_bigt_bigmat, 3},
    {"_simer_impute_marker", (DL_FUNC) &_simer_impute_marker, 4},
    {"_simer_hasNA", (DL_FUNC) &_simer_hasNA, 5},
    {"_simer_hasNABed", (DL_FUNC) &_simer_hasNABed, 5},
    {"_simer_PedigreeCorrector", (DL_FUNC) &_simer_PedigreeCorrector, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_simer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
