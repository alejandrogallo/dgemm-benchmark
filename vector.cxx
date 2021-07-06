#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <chrono>
#include <memory>
#include "utils.h"
#include "dgemm.h"

int main(int argc, char ** argv){

  int rank, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  const bool vector(hauta::option<bool>(argc, argv, "--vector", false))
           , initialize(hauta::option<bool>(argc, argv, "--init", false))
           ;
  const std::string container = vector ? "std::vector" : "double[]";

  const
  size_t No = hauta::option<size_t>(argc, argv, "--no")
       , Nv = hauta::option<size_t>(argc, argv, "--nv")
       , iterations = hauta::option<size_t>(argc, argv, "-i", 1)
       , avg = hauta::option<size_t>(argc, argv, "--avg", 1)
       , m = No
       , n = No*No
       , k = Nv
       ;

  Timings chrono;
  Averages averages;

  double *A, *B, *C;

  //std::vector<double> vA, vB, vC;
  //vA.reserve(m*k);
  //vB.reserve(n*k);
  //vC.reserve(n*m);

  //std::vector<double> vA(m*k), vB(n*k), vC(n*m);
  std::vector<double> vA, vB, vC;

  vA.reserve(m*k);
  vB.reserve(n*k);
  vC.reserve(n*m);

  if (initialize) {
    vA.resize(m*k);
    vB.resize(n*k);
    vC.resize(n*m);
  }

  for (auto &i: vA) i += 5;
  for (auto &i: vB) i += 5;
  for (auto &i: vC) i += 5;


  if (!vector) {
    A = new double[m*k];
    B = new double[k*n];
    C = new double[m*n];
    if (initialize) {
      for (size_t i=0; i<m*k; i++) A[i] = 0.0;
      for (size_t i=0; i<k*n; i++) B[i] = 0.0;
      for (size_t i=0; i<m*n; i++) C[i] = 0.0;
    }
  } else {
    A = vA.data();
    B = vB.data();
    C = vC.data();
  }

  const double flopCount
    = double(m * n * k)
    * 6
    * 2
    * double(iterations)
    / 1e9
    ;

  LOG << "======= VECTOR ======\n";
  LOG << SHOW_VAR(container) << "\n";
  LOG << SHOW_VAR(initialize) << "\n";
  LOG << SHOW_VAR(np) << "\n";
  LOG << SHOW_VAR(No) << "\n";
  LOG << SHOW_VAR(Nv) << "\n";
  LOG << SHOW_VAR(iterations) << "\n";
  LOG << SHOW_VAR(flopCount) << "\n";
  LOG << SHOW_MACRO(GIT_COMMIT) << "\n";
  LOG << SHOW_MACRO(DATE) << "\n";
  LOG << SHOW_MACRO(CONFIG) << "\n";
  LOG << SHOW_MACRO(COMPILER_VERSION) << "\n";
  LOG << SHOW_VAR(m) << " " << SHOW_VAR(n) << " " << SHOW_VAR(k) << "\n";


  for (size_t __avg = 1; __avg <= avg ; __avg++) {

    chrono["main"].start();
    for (size_t it = 0; it < iterations; it++) {
      for (size_t _m = 0; _m < m; _m++)
      for (size_t _n = 0; _n < n; _n++)
      for (size_t _k = 0; _k < k; _k++) {
        C[_n + _m * n] = A[_k + _m * k] * B[_n + _k * n] + C[_n + _m * n];
      }
    }
    chrono["main"].stop();
    averages["flops:main"].push(__avg * flopCount / chrono["main"].count());

  }

  LOG << "main: " << chrono["main"].count() << std::endl;

  for (auto const& a: averages)
    LOG << a.first << " "
        << a.second.count() << " ± " << a.second.sigma()
        << " :avg " << a.second.size()
        << "\n"
        ;

  MPI_Finalize();
  return 0;

}

