#include <iostream>
#include <vector>
#include <iomanip>
#include <mpi.h>
#include "utils.h"
#include "dgemm.h"


int main(int argc, char ** argv){

  int rank, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Timings chrono;
  double one(1.0), m_one(-1.0), zero(0.0);

  size_t No = hauta::option<size_t>(argc, argv, "--no")
       , Nv = hauta::option<size_t>(argc, argv, "--nv")
       , iterations = hauta::option<size_t>(argc, argv, "-i", 1)
       ;

  const int NoNo = No*No, NoNoNo = No*No*No;
  const double flopCountParticles = double(NoNoNo) * double(Nv) * 2.0 * 6.0 / 1e9;
  const double flopCountHoles = double(NoNoNo) * double(No) * 2.0 * 6.0 / 1e9;
  const double flopCount = flopCountHoles + flopCountParticles;

  std::vector<double>
      TABhh, VhhhC, Tijk
    , TAphh, VBCph
    ;

  TABhh.reserve(No*No);
  VhhhC.reserve(NoNoNo);
  Tijk.reserve(NoNoNo);
  TAphh.reserve(NoNo*Nv);
  VBCph.reserve(Nv*No);

  LOG << "∷∷∷∷∷∷∷ ĐG∃MM ======\n";
  LOG << SHOW_VAR(np) << "\n";
  LOG << SHOW_VAR(No) << "\n";
  LOG << SHOW_VAR(Nv) << "\n";
  LOG << SHOW_VAR(iterations) << "\n";
  LOG << SHOW_VAR(flopCount) << "\n";
  LOG << SHOW_VAR(flopCountHoles) << "\n";
  LOG << SHOW_VAR(flopCountParticles) << "\n";
#if defined(BLIS_ARCH)
  LOG << SHOW_MACRO(BLIS_ARCH) << "\n";
#endif
  LOG << SHOW_MACRO(GIT_COMMIT) << "\n";
  LOG << SHOW_MACRO(DATE) << "\n";
  LOG << SHOW_MACRO(CONFIG) << "\n";
  LOG << SHOW_MACRO(COMPILER_VERSION) << "\n";

  chrono["doubles"].start();
  for (size_t it = 0; it < iterations; it++) {
    chrono["start:stop"].start();
    chrono["start:stop"].stop();

    // holes
    chrono["holes"].start();
    for (size_t i = 0; i < 6; i++) {
      chrono["holes:dgemm"].start();
      dgemm_( "N"
            , "T"
            , (int const*)&NoNo
            , (int const*)&No
            , (int const*)&No
            , &m_one
            , VhhhC.data()
            , (int const*)&NoNo
            , TABhh.data()
            , (int const*)&No
            , &zero
            , Tijk.data()
            , (int const*)&NoNo
            );
      chrono["holes:dgemm"].stop();
#if !defined(NO_REORDER)
      chrono["holes:reorder"].start();
      for (int ijk = 0; ijk < NoNoNo; ijk++) {
        Tijk.data()[ijk] = 1.0;
      }
      chrono["holes:reorder"].stop();
#endif
    }
    chrono["holes"].stop();

    // particles
    chrono["particles"].start();
    for (size_t i = 0; i < 6; i++) {
      chrono["particles:dgemm"].start();
      dgemm_( "T"
            , "N"
            , (int const*)&NoNo
            , (int const*)&No
            , (int const*)&Nv
            , &one
            , TAphh.data()
            , (int const*)&Nv
            , VBCph.data()
            , (int const*)&Nv
            , &zero
            , Tijk.data()
            , (int const*)&NoNo
            );
      chrono["particles:dgemm"].stop();
#if !defined(NO_REORDER)
      chrono["particles:reorder"].start();
      for (int ijk = 0; ijk < NoNoNo; ijk++) {
        Tijk.data()[ijk] = 1.0;
      }
      chrono["particles:reorder"].stop();
#endif
    }
    chrono["particles"].stop();

  }
  chrono["doubles"].stop();




  // PRINT TIMINGS
  for (auto const& pair: chrono)
    LOG << std::setprecision(6) << std::setw(6)
        << pair.first << ": "
        << pair.second.count()
        << std::endl;

  LOG
    << "flops:doubles "
    << flopCount * iterations / chrono["doubles"].count()
    << "\n";
  LOG
    << "flops:doubles:no-reorder "
    << flopCount * iterations
    / (chrono["particles:dgemm"].count() + chrono["holes:dgemm"].count())
    << "\n";
  LOG
    << "flops:particles "
    << flopCountParticles * iterations / chrono["particles:dgemm"].count()
    << "\n";
  LOG
    << "flops:holes "
    << flopCountHoles * iterations / chrono["holes:dgemm"].count()
    << "\n";

  MPI_Finalize();

  return 0;

}

