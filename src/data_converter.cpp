#include "simer_omp.h"
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <progress.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;

// ***** BFILE *****

template <typename T>
void write_bfile(XPtr<BigMatrix> pMat, std::string bed_file, double NA_C, int threads=0, bool verbose=true) {
  // check input
  string ending = ".bed";
  if (bed_file.length() <= ending.length() ||
      0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending)) {
    bed_file += ending;
  }
  
  // define
  T c;
  omp_setup(threads);
  int m = pMat->nrow();
  int nind = pMat->ncol();
  int n = pMat->ncol() / 4;  // 4 individual = 1 bit
  if (pMat->ncol() % 4 != 0) 
    n++;
  
  vector<uint8_t> geno(n);
  MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
  FILE *fout;
  fout = fopen(bed_file.c_str(), "wb");
  
  // progress bar
  Progress progress(m, verbose);
  
  // magic number of bfile
  const unsigned char magic_bytes[] = { 0x6c, 0x1b, 0x01 };
  fwrite((char*)magic_bytes, 1, 3, fout);
  
  // map
  std::map<T, int> code;
  code[0] = 3;
  code[1] = 2;
  code[2] = 0;
  code[static_cast<T>(NA_C)] = 1;
  
  // write bfile
  for (int i = 0; i < m; i++) {
  #pragma omp parallel for private(c)
    for (int j = 0; j < n; j++) {
      uint8_t p = 0;
      for (int x = 0; x < 4 && (4 * j + x) < nind; x++) {
        c = mat[4 * j + x][i];
        p |= code[c] << (x*2);
      }
      geno[j] = p;
    }
    fwrite((char*)geno.data(), 1, geno.size(), fout);
    progress.increment();
  }
  fclose(fout);
  return;
}

// [[Rcpp::export]]
void write_bfile(SEXP pBigMat, std::string bed_file, int threads=0, bool verbose=true) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return write_bfile<char>(xpMat, bed_file, NA_CHAR, threads, verbose);
  case 2:
    return write_bfile<short>(xpMat, bed_file, NA_SHORT, threads, verbose);
  case 4:
    return write_bfile<int>(xpMat, bed_file, NA_INTEGER, threads, verbose);
  case 8:
    return write_bfile<double>(xpMat, bed_file, NA_REAL, threads, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

template <typename T>
void read_bfile(std::string bed_file, XPtr<BigMatrix> pMat, long maxLine, double NA_C, int threads=0, bool verbose=true) {
  // check input
  if (!boost::ends_with(bed_file, ".bed")) {
    bed_file += ".bed";
  }
  
  // define
  omp_setup(threads);
  size_t ind = pMat->ncol();
  long n = ind / 4;  // 4 individual = 1 bit
  if (ind % 4 != 0) 
    n++; 
  char *buffer;
  long buffer_size;
  MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
  
  // map
  std::map<int, T> code;
  code[3] = 0;
  code[2] = 1;
  code[1] = static_cast<T>(NA_C);
  code[0] = 2;
  
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
  Progress progress(n_block, verbose);
  
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
      
      for (size_t x = 0; x < 4 && (c + x) < ind; x++) {
        mat[c + x][r] = code[(p >> (2*x)) & 0x03];
      }
    }
    progress.increment();
  }
  fclose(fin);
  return;
}

// [[Rcpp::export]]
void read_bfile(std::string bed_file, SEXP pBigMat, long maxLine, int threads=0, bool verbose=true) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return read_bfile<char>(bed_file, xpMat, maxLine, NA_CHAR, threads, verbose);
  case 2:
    return read_bfile<short>(bed_file, xpMat, maxLine, NA_SHORT, threads, verbose);
  case 4:
    return read_bfile<int>(bed_file, xpMat, maxLine, NA_INTEGER, threads, verbose);
  case 8:
    return read_bfile<double>(bed_file, xpMat, maxLine, NA_REAL, threads, verbose);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}
