#include "simer_omp.h"
#include <boost/algorithm/string.hpp>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigmemory, BH)]]
using namespace std;
using namespace Rcpp;

template <typename T>
bool hasNA(XPtr<BigMatrix> pMat, double NA_C, const int threads=0) {
  
  omp_setup(threads);
  size_t m = pMat->nrow();
  size_t n = pMat->ncol();
  bool HasNA = false;
  
  MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
  #pragma omp parallel for schedule(dynamic) shared(HasNA)
  for (size_t j = 0; j < n; j++) {
    if(HasNA)   continue;
    for (size_t i = 0; i < m; i++) {
      if (mat[j][i] == NA_C) {
        HasNA = true;
      }
    }
  }
  return HasNA;
}

// [[Rcpp::export]]
bool hasNA(SEXP pBigMat, const int threads=0) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return hasNA<char>(xpMat, NA_CHAR, threads);
  case 2:
    return hasNA<short>(xpMat, NA_SHORT, threads);
  case 4:
    return hasNA<int>(xpMat, NA_INTEGER, threads);
  case 8:
    return hasNA<double>(xpMat, NA_REAL, threads);
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
  long n = ind / 4;  // 4 individual = 1 bit
  if (ind % 4 != 0) 
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
  
  // loop file
  size_t cond;
  long block_start;
  for (int i = 0; i < n_block; i++) {
    buffer = new char [buffer_size];
    n_bytes_read = static_cast<size_t>(fread(buffer, 1, buffer_size, fin));
    
    // i: current block, j: current bit.
    block_start = i * buffer_size;
    cond = min(buffer_size, length - 3 - block_start);
    #pragma omp parallel for schedule(static)
    for (size_t j = 0; j < cond; j++) {
      // bit -> item in matrix
      size_t r = (block_start + j) / n;
      size_t c = (block_start + j) % n * 4;
      uint8_t p = buffer[j];
      if (HasNA) continue;
      for (size_t x = 0; x < 4 && (c + x) < ind; x++) {
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
