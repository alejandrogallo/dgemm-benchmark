#include <iostream>
#include <vector>
#include <iomanip>
#include <mpi.h>
#include "utils.h"
#include "dgemm.h"


int main(int argc, char ** argv){
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Timings chrono;
  double one(1.0), m_one(-1.0), zero(0.0);
  const int
      No(atoi(argv[1]))
    , Nv(atoi(argv[2]))
    , iterations(atoi(argv[3]))
    ;

  const int NoNo = No*No, NoNoNo = No*No*No;
  const double flopCountParticles = double(NoNoNo) * double(Nv) * 2.0 * 6.0 / 1e9;
  const double flopCountHoles = double(NoNoNo) * double(No) * 2.0 * 6.0 / 1e9;
  const double flopCount = flopCountHoles + flopCountParticles;

  std::vector<double>
      TABhh(No*No), VhhhC(NoNoNo), Tijk(NoNoNo)
    , TAphh(NoNo*Nv), VBCph(Nv*No);

  LOG << "Doing DGEMM Tests\n";
  LOG << "»»»»»»»»»»»»»»»»»\n";
  LOG << SHOW_VAR(No) << "\n";
  LOG << SHOW_VAR(Nv) << "\n";
  LOG << SHOW_VAR(iterations) << "\n";
  LOG << SHOW_VAR(flopCount) << "\n";
  LOG << SHOW_VAR(flopCountHoles) << "\n";
  LOG << SHOW_VAR(flopCountParticles) << "\n";
#if defined(HAS_INTEL)
  LOG << "intel compiler\n";
#elif defined(HAS_GCC)
  LOG << "gcc compiler\n";
#endif
#if defined(BLIS_ARCH)
  LOG << SHOW_MACRO(BLIS_ARCH) << "\n";
#endif
  LOG << SHOW_MACRO(GIT_COMMIT) << "\n";
  LOG << SHOW_MACRO(CONFIG) << "\n";
  LOG << SHOW_MACRO(COMPILER_VERSION) << "\n";

  chrono["doubles"].start();
  for (int it = 0; it < iterations; it++) {
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

