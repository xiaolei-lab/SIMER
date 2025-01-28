// Data pre-processing module
// 
// Copyright (C) 2016-2018 by Xiaolei Lab
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

// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

// ***** BFILE *****

template <typename T>
void write_bfile(XPtr<BigMatrix> pMat, std::string bed_file, double NA_C, bool mrkbycol = true, int threads=0, bool verbose=true) {
    // check input
    string ending = ".bed";
    if (bed_file.length() <= ending.length() ||
        0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending)) {
        bed_file += ending;
    }
    
    // define
    T c;
    omp_setup(threads);
    int m = (mrkbycol ? pMat->ncol() : pMat->nrow());
    int nind = (mrkbycol ? pMat->nrow() : pMat->ncol());
    int n = nind / 4;  // 4 individual = 1 bit
    if (nind % 4 != 0) 
        n++;
    
    vector<uint8_t> geno(n);
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    FILE *fout;
    fout = fopen(bed_file.c_str(), "wb");
    
    // progress bar
    MinimalProgressBar_perc pb;
    Progress progress(m, verbose, pb);
    
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
    if (mrkbycol) {
        for (int i = 0; i < m; i++) {
            #pragma omp parallel for private(c)
            for (int j = 0; j < n; j++) {
                uint8_t p = 0;
                for (int x = 0; x < 4 && (4 * j + x) < nind; x++) {
                    c = mat[i][4 * j + x];
                    p |= code[c] << (x*2);
                }
                geno[j] = p;
            }
            fwrite((char*)geno.data(), 1, geno.size(), fout);
            progress.increment();
        }
    }else{
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
    }
    fclose(fout);
    return;
}

// [[Rcpp::export]]
void write_bfile(SEXP pBigMat, std::string bed_file, bool mrkbycol = true, int threads=0, bool verbose=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return write_bfile<char>(xpMat, bed_file, NA_CHAR, mrkbycol, threads, verbose);
    case 2:
        return write_bfile<short>(xpMat, bed_file, NA_SHORT, mrkbycol, threads, verbose);
    case 4:
        return write_bfile<int>(xpMat, bed_file, NA_INTEGER, mrkbycol, threads, verbose);
    case 8:
        return write_bfile<double>(xpMat, bed_file, NA_REAL, mrkbycol, threads, verbose);
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
    size_t ind = pMat->nrow();
    long n = ind / 4;  // 4 individual = 1 bit
    if (ind % 4 != 0) 
        n++; 
    char *buffer;
    long buffer_size;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    // map
    std::map<int, T> code;
    code[3] = static_cast<T>(0);
    code[2] = static_cast<T>(1);
    code[1] = static_cast<T>(NA_C);
    code[0] = static_cast<T>(2);
    
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
    MinimalProgressBar_perc pb;
    Progress progress(n_block, verbose, pb);
    
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
            size_t r = j / n + i * maxLine;
            size_t c = j % n * 4;
            uint8_t p = buffer[j];
            
            for (size_t x = 0; x < 4 && (c + x) < ind; x++) {
                mat[r][c + x] = code[(p >> (2*x)) & 0x03];
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