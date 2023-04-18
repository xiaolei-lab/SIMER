#ifndef MINIMALPROGRESSBAR_H_
#define MINIMALPROGRESSBAR_H_

#include <Rcpp.h>
#include <R_ext/Print.h>
#include <progress.hpp>
#include "progress_bar.hpp"

// [[Rcpp::plugins(cpp11)]]
using namespace std;
using namespace Rcpp;

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

#endif
