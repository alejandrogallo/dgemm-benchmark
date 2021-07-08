#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

extern
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

typedef enum Initialization { NONE = 0, LOOP = 7, CALLOC = 5 }
  Initialization;

Initialization initialization_from_env() {
  if (getenv("INIT") != NULL && strcmp(getenv("INIT"), "LOOP") == 0) {
    return LOOP;
  } else if (getenv("INIT") != NULL && strcmp(getenv("INIT"), "CALLOC") == 0) {
    return CALLOC;
  } else {
    return NONE;
  }
}

int main(int argc, char ** argv){

  MPI_Init(0, 0);
  if (argc < 4) {
    printf("n, m, k - matrix dimensions\n");
    return 1;
  }

  const Initialization initialize = initialization_from_env();
  const int use_dgemm = strcmp(getenv("USE_DGEMM"), "JA") == 0;

  const int m = atoi(argv[1]);
  const int n = atoi(argv[2]);
  const int k = atoi(argv[3]);
  int _m, _n, _k, i;

  double *A, *B, *C;
  double one = (1.0), tmp;

  if (initialize == CALLOC) {
    A = (double*)calloc(m*k, sizeof(double));
    B = (double*)calloc(k*n, sizeof(double));
    C = (double*)calloc(m*n, sizeof(double));
    printf("calloc flops\n");
  } else {
    A = (double*)malloc(sizeof(double)*m*k);
    B = (double*)malloc(sizeof(double)*k*n);
    C = (double*)malloc(sizeof(double)*m*n);
  }

  if (initialize == LOOP) {
    printf("initializing \n");
    for (_m=0; _m<m*k; _m++) A[_m] = 897.0;
    for (_m=0; _m<k*n; _m++) B[_m] = 897.0;
    for (_m=0; _m<m*n; _m++) C[_m] = 89.0;
  }

  const double flops = (2*k*n*m);
  dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);

  double start = MPI_Wtime();

  if (use_dgemm) {
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
  } else {
#define MYDGEMM                  \
    for (_m = 0; _m < m; _m++) { \
    for (_n = 0; _n < n; _n++) { \
      tmp = C[_n + _m * n];      \
      C[_n + _m * n] = 46.0;     \
      A[_m] = 46.0;              \
      B[_n] = 46.0;              \
    for (_k = 0; _k < k; _k++) { \
      tmp += A[_k + _m * k]      \
           * B[_k + _n * k]      \
           ;                     \
    }                            \
      C[_n + _m * n] = tmp;      \
    }}

    for (i = 0; i < 10; i++) {
      MYDGEMM
      MYDGEMM
      MYDGEMM
    }
#undef MYDGEMM
  }

  double end = MPI_Wtime();
  double  est = end - start;
  printf("time: %f->%f : %f\n", start, end, est);
  printf("m %d n %d k %d (%d) flops: %.1e %f\n",
         m, n, k,
         initialize,
         flops,
         (use_dgemm ? 3 : 30)
           * flops
           / 1e9
           / est);

  MPI_Finalize();
  return 0;

}

