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

  Timings chrono;
  auto chrono_main = chrono["main"];
  int m(hauta::option<int>(argc, argv, "-m"))
    , n(hauta::option<int>(argc, argv, "-n"))
    , k(hauta::option<int>(argc, argv, "-k"))
    , iterations(hauta::option<int>(argc, argv, "-i", 1))
    ;
  const bool warmup(hauta::option<bool>(argc, argv, "-/warmup", true))
           , vector(hauta::option<bool>(argc, argv, "--vector", false))
           ;

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
  std::chrono::duration<double> duration;

  LOG << SHOW_VAR(flopCount) << "\n";
  LOG << SHOW_VAR(iterations) << "\n";
  LOG << SHOW_VAR(np) << "\n";

  if (warmup) {
    LOG << "Warming up \n";
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t it = 0; it < iterations; it++) {
    chrono_main.start();
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    dgemm_("N", "N", &m, &n, &k, &one, A, &m, B, &k, &one, C, &m);
    chrono_main.stop();
  }
  auto end = std::chrono::high_resolution_clock::now();

  double est = 0.001 *std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  duration = end - start;
  LOG << "main: " << chrono_main.count() << std::endl;
  LOG << "duration: " << duration.count() << std::endl;
  LOG << "time: " << est << std::endl;
  LOG << m << " "
      << n << " "
      << k << " flops: "
      << flopCount << " "<< 3*flopCount/1e9/est << std::endl;
  LOG << m << " "
      << n << " "
      << k
      << " flops: " << flopCount << " "<< 3*flopCount/1e9/(chrono_main.count()) << std::endl;

  MPI_Finalize();
  return 0;

}

