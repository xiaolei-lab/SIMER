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
  MinimalProgressBar(const string str = "Calculating in process")  {
    _finalized = false;
    _str = str;
  }
  ~MinimalProgressBar() {}
  void display() {}
  void update(float progress) {
    if (_finalized) return;
    REprintf("\r");
    REprintf(_str.c_str());
    REprintf("...finished %u%%", (int)(progress * 100));
  }
  void end_display() {
    if (_finalized) return;
    REprintf("\r");
    REprintf(_str.c_str());
    REprintf("...[finished 100%%]");
    REprintf("\n");
    _finalized = true;
  }
  private:
  bool _finalized;
  string _str;
};

#endif
