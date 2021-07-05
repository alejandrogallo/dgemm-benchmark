#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <chrono>
#include "utils.h"
#include "dgemm.h"

int main(int argc, char ** argv){

  int rank, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  const bool warmup(hauta::option<bool>(argc, argv, "-/warmup", true))
           , vector(hauta::option<bool>(argc, argv, "--vector", false))
           , holes(hauta::option<bool>(argc, argv, "--holes", false))
           ;

  const
  size_t No = hauta::option<size_t>(argc, argv, "--no")
       , Nv = hauta::option<size_t>(argc, argv, "--nv")
       , iterations = hauta::option<size_t>(argc, argv, "-i", 1)
       ;

  int m = holes ? No*No : No
    , n = holes ? No    : No*No
    , k = holes ? No    : Nv
    ;

  Timings chrono;
  auto chrono_main = chrono["main"];

  double *A, *B, *C;
  std::vector<double> vA, vB, vC;
  vA.reserve(m*k);
  vB.reserve(n*k);
  vC.reserve(n*m);
  if (!vector) {
    LOG << "Using carrays\n";
    A = new double[m*k];
    B = new double[k*n];
    C = new double[m*n];
  } else {
    LOG << "Using std::vector\n";
    A = vA.data();
    B = vB.data();
    C = vC.data();
  }

  double one(1.0);
  const double flopCount(double(n*m*2*k) * double(iterations));

  LOG << "======= BLAS ======\n";
  LOG << SHOW_VAR(np) << "\n";
  LOG << SHOW_VAR(No) << "\n";
  LOG << SHOW_VAR(Nv) << "\n";
  LOG << SHOW_VAR(iterations) << "\n";
  LOG << SHOW_VAR(flopCount) << "\n";
  LOG << SHOW_VAR(holes) << "\n";
#if defined(BLIS_ARCH)
  LOG << SHOW_MACRO(BLIS_ARCH) << "\n";
#endif
  LOG << SHOW_MACRO(GIT_COMMIT) << "\n";
  LOG << SHOW_MACRO(DATE) << "\n";
  LOG << SHOW_MACRO(CONFIG) << "\n";
  LOG << SHOW_MACRO(COMPILER_VERSION) << "\n";
  LOG << SHOW_VAR(m) << " " << SHOW_VAR(n) << " " << SHOW_VAR(k) << "\n";


  if (warmup) {
    LOG << "Warming up \n";
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
  }

  chrono_main.start();
  for (size_t it = 0; it < iterations; it++) {
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
  }
  chrono_main.stop();

  LOG << "main: " << chrono_main.count() << std::endl;
  LOG << "flops: " << 6 * flopCount / 1e9 / chrono_main.count()
      << std::endl;

  MPI_Finalize();
  return 0;

}

