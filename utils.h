#pragma once
#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <map>
#include <string>
#define Q(x) #x
#define QUOTE(x) Q(x)
#define SHOW_MACRO(x) #x " = " QUOTE(x)
#define SHOW_VAR(x) #x " = " << x

#if defined(HAS_MKL) || defined(HAS_LAPACK)

extern "C" {
  void dgemm_(
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    double *alpha,
    const double *A,
    const int *lda,
    const double *B,
    const int *ldb,
    double *beta,
    double *C,
    const int *ldc
  );
}

#endif

struct Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Event = std::chrono::time_point<Clock>;
  std::chrono::duration<double> duration;
  Event _start;
  Event start() { return _start = Clock::now(); }
  Event stop() {
    Event const _end = Clock::now();
    duration += _end - _start;
    return _end;
  }
  double count() const { return duration.count(); }
};

using Timings = std::map<std::string, Timer>;
