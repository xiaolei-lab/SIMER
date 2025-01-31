// Copyright (C) 2016-2019 by Xiaolei Team
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "simer.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

template <typename T>
void impute_marker(XPtr<BigMatrix> pMat, bool mrkbycol = true, int threads=0, bool verbose=true) {
    omp_setup(threads);

    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    const size_t n = (mrkbycol ? pMat->nrow() : pMat->ncol());
    const size_t m = (mrkbycol ? pMat->ncol() : pMat->nrow());
    
    MinimalProgressBar_perc pb;
    Progress progress(m, verbose, pb);
    
    // loop marker
    if(mrkbycol){
        #pragma omp parallel for
        for (size_t i = 0; i < m; i++) {
            std::vector<size_t> na_index = {};;
            size_t counts[3] = { 0 };
            
            // count allele, record missing index 
            for (size_t j = 0; j < n; j++) {
                switch(int(mat[i][j])) {
                case 0: counts[0]++; break;
                case 1: counts[1]++; break;
                case 2: counts[2]++; break;
                default: na_index.push_back(j);
                }
            }
            
            // find major allele
            T major = counts[2] > counts[1] ? (counts[2] > counts[0] ? 2 : 0) : (counts[1] > counts[0] ? 1 : 0);
            
            // impute
            for (auto&& x: na_index) {
                mat[i][x] = major;   
            }
            progress.increment();
        }
    }else{
        #pragma omp parallel for
        for (size_t i = 0; i < m; i++) {
            std::vector<size_t> na_index = {};;
            size_t counts[3] = { 0 };
            
            // count allele, record missing index 
            for (size_t j = 0; j < n; j++) {
                switch(int(mat[j][i])) {
                case 0: counts[0]++; break;
                case 1: counts[1]++; break;
                case 2: counts[2]++; break;
                default: na_index.push_back(j);
                }
            }
            
            // find major allele
            T major = counts[2] > counts[1] ? (counts[2] > counts[0] ? 2 : 0) : (counts[1] > counts[0] ? 1 : 0);
            
            // impute
            for (auto&& x: na_index) {
                mat[x][i] = major;   
            }
            progress.increment();
        }
    }
}

// [[Rcpp::export]]
void impute_marker(SEXP pBigMat, bool mrkbycol = true, int threads=0, bool verbose=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return impute_marker<char>(xpMat, mrkbycol, threads, verbose);
    case 2:
        return impute_marker<short>(xpMat, mrkbycol, threads, verbose);
    case 4:
        return impute_marker<int>(xpMat, mrkbycol, threads, verbose);
    case 8:
        return impute_marker<double>(xpMat, mrkbycol, threads, verbose);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

template <typename T>
bool hasNA(XPtr<BigMatrix> pMat, bool mrkbycol = true, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const int threads = 1) {
    
    omp_setup(threads);

    bool HasNA = false;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    if(geno_ind.isNotNull()){
        uvec _geno_ind = as<uvec>(geno_ind) - 1;
        int n = _geno_ind.n_elem;
        if(marker_ind.isNotNull()){
            uvec _marker_ind = as<uvec>(marker_ind) - 1;
            int m = _marker_ind.n_elem;
            if(mrkbycol){
                #pragma omp parallel for shared(HasNA)
                for (int j = 0; j < m; j++) {
                    if(HasNA)   continue;
                    for (int i = 0; i < n; i++) {
                        // if (mat[_marker_ind[j]][_geno_ind[i]] == NA_C) {
                        if (isna(mat[_marker_ind[j]][_geno_ind[i]])) {
                            HasNA = true;
                        }
                    }
                }
            }else{
                #pragma omp parallel for shared(HasNA)
                for (int j = 0; j < n; j++) {
                    if(HasNA)   continue;
                    for (int i = 0; i < m; i++) {
                        // if (mat[_geno_ind[j]][_marker_ind[i]] == NA_C) {
                        if (isna(mat[_geno_ind[j]][_marker_ind[i]])) {
                            HasNA = true;
                        }
                    }
                }
            }
        }else{
            if(mrkbycol){
                #pragma omp parallel for shared(HasNA)
                for (int j = 0; j < pMat->ncol(); j++) {
                    if(HasNA)   continue;
                    for (int i = 0; i < n; i++) {
                        // if (mat[j][_geno_ind[i]] == NA_C) {
                        if (isna(mat[j][_geno_ind[i]])) {
                            HasNA = true;
                        }
                    }
                }
            }else{
                #pragma omp parallel for shared(HasNA)
                for (int j = 0; j < n; j++) {
                    if(HasNA)   continue;
                    for (int i = 0; i < pMat->nrow(); i++) {
                        // if (mat[_geno_ind[j]][i] == NA_C) {
                        if (isna(mat[_geno_ind[j]][i])) {
                            HasNA = true;
                        }
                    }
                }
            }
        }
    }else{
        if(marker_ind.isNotNull()){
            uvec _marker_ind = as<uvec>(marker_ind) - 1;
            int m = _marker_ind.n_elem;
            if(mrkbycol){
                #pragma omp parallel for shared(HasNA)
                for (int j = 0; j < m; j++) {
                    if(HasNA)   continue;
                    for (int i = 0; i < pMat->nrow(); i++) {
                        // if (mat[_marker_ind[j]][i] == NA_C) {
                        if (isna(mat[_marker_ind[j]][i])) {
                            HasNA = true;
                        }
                    }
                }
            }else{
                #pragma omp parallel for shared(HasNA)
                for (int j = 0; j < pMat->ncol(); j++) {
                    if(HasNA)   continue;
                    for (int i = 0; i < m; i++) {
                        // if (mat[j][_marker_ind[i]] == NA_C) {
                        if (isna(mat[j][_marker_ind[i]])) {
                            HasNA = true;
                        }
                    }
                }
            }
        }else{
            #pragma omp parallel for shared(HasNA)
            for (int j = 0; j < pMat->ncol(); j++) {
                if(HasNA)   continue;
                for (int i = 0; i < pMat->nrow(); i++) {
                    // if (mat[j][i] == NA_C) {
                    if (isna(mat[j][i])) {
                        HasNA = true;
                    }
                }
            }
        }
    }
    return HasNA;
}

// [[Rcpp::export]]
bool hasNA(SEXP pBigMat, bool mrkbycol = true, const Nullable<arma::uvec> geno_ind = R_NilValue, const Nullable<arma::uvec> marker_ind = R_NilValue, const int threads=1) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return hasNA<char>(xpMat, mrkbycol, geno_ind, marker_ind, threads);
    case 2:
        return hasNA<short>(xpMat, mrkbycol, geno_ind, marker_ind, threads);
    case 4:
        return hasNA<int>(xpMat, mrkbycol, geno_ind, marker_ind, threads);
    case 8:
        return hasNA<double>(xpMat, mrkbycol, geno_ind, marker_ind, threads);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

// [[Rcpp::export]]
bool hasNABed(std::string bed_file, int ind, long maxLine, int threads=0, bool verbose=true) {
  // check input
  if (!boost::ends_with(bed_file, ".bed")) {
    bed_file += ".bed";
  }
  
  // define
  omp_setup(threads);
  size_t indt = ind;
  long n = indt / 4;  // 4 individual = 1 bit
  if (indt % 4 != 0) 
    n++; 
  char *buffer;
  long buffer_size;
  bool HasNA = false;
  
  // open file
  FILE *fin;
  fin = fopen(bed_file.c_str(), "rb");
  fseek(fin, 0, SEEK_END);
  long length = ftell(fin);
  rewind(fin);
  
  // get buffer_size
  buffer_size = maxLine > 0 ? (maxLine * n) : (length - 3);
  
  // progress bar
  int n_block = (length - 3) / buffer_size;
  if ((length - 3) % buffer_size != 0) { n_block++; }
  
  // magic number of bfile
  buffer = new char [3];
  size_t n_bytes_read = static_cast<size_t>(fread(buffer, 1, 3, fin));
  if (n_bytes_read != 3) {    // new if
    Rcpp::stop("It is not a normal binary file!");
  }

  // loop file
  size_t cond;
  long block_start;
  for (int i = 0; i < n_block; i++) {
    buffer = new char [buffer_size];
    n_bytes_read = static_cast<size_t>(fread(buffer, 1, buffer_size, fin));
    
    // i: current block, j: current bit.
    block_start = i * buffer_size;
    cond = min(buffer_size, length - 3 - block_start);

    #pragma omp parallel for
    for (size_t j = 0; j < cond; j++) {
      // bit -> item in matrix
      // size_t r = j / n + i * maxLine;
      size_t c = j % n * 4;
      uint8_t p = buffer[j];

      if (HasNA) continue;

      for (size_t x = 0; x < 4 && (c + x) < indt; x++) {
        if (1 == ((p >> (2*x)) & 0x03)) {
          HasNA = true;
          break;
        }
      }
    }
  }
  fclose(fin);
  
  return HasNA;
}
