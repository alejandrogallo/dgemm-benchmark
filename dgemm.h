#pragma once

#ifdef HAS_BLIS
#include <blis.h>
#endif

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
