#include <iostream>
#include <vector>
#include <iomanip>
#include "utils.h"
#include "dgemm.h"
#include <mpi.h>

#define _IJK_(i, j, k) i + j*No + k*NoNo
#define REORDER(__II, __JJ, __KK)                                 \
  chrono["doubles:reorder"].start();                              \
  for (size_t k = 0; k < No; k++)                                 \
  for (size_t j = 0; j < No; j++)                                 \
  for (size_t i = 0; i < No; i++) {                               \
    Tijk[_IJK_(i, j, k)] += _t_buffer[_IJK_(__II, __JJ, __KK)];   \
  }                                                               \
  chrono["doubles:reorder"].stop();
#define DGEMM_PARTICLES(__A, __B)    \
         dgemm_( "T"                 \
               , "N"                 \
               , (int const*)&NoNo   \
               , (int const*)&No     \
               , (int const*)&Nv     \
               , &one                \
               , __A                 \
               , (int const*)&Nv     \
               , __B                 \
               , (int const*)&Nv     \
               , &zero               \
               , _t_buffer.data()    \
               , (int const*)&NoNo   \
               );
#define DGEMM_HOLES(__A, __B, __TRANSB)  \
         dgemm_( "N"                     \
               , __TRANSB                \
               , (int const*)&NoNo       \
               , (int const*)&No         \
               , (int const*)&No         \
               , &m_one                  \
               , __A                     \
               , (int const*)&NoNo       \
               , __B                     \
               , (int const*)&No         \
               , &zero                   \
               , _t_buffer.data()        \
               , (int const*)&NoNo       \
               );

void doDoubles(size_t No, size_t Nv,
               double const* VhhhC, double const* TABhh,
               double const* TAphh, double const* VBCph,
               double* Tijk,
               Timings& chrono) {

  const int NoNo = No*No, NoNoNo = No*No*No;
  double one(1.0), m_one(-1.0), zero(0.0);
  std::vector<double> _t_buffer(NoNoNo);

  chrono["doubles:holes"].start();
  { // Holes part ============================================================
    // VhhhC[i + k*No + L*NoNo] * TABhh[L + j*No]; H1
    DGEMM_HOLES(VhhhC, TABhh, "N")
    REORDER(i, k, j)
    // VhhhC[j + k*No + L*NoNo] * TABhh[i + L*No]; H0
    DGEMM_HOLES(VhhhC, TABhh, "T")
    REORDER(j, k, i)
    // VhhhB[i + j*No + L*NoNo] * TAChh[L + k*No]; H5
    DGEMM_HOLES(VhhhC, TABhh, "N")
    REORDER(i, j, k)
    // VhhhB[k + j*No + L*NoNo] * TAChh[i + L*No]; H3
    DGEMM_HOLES(VhhhC, TABhh, "T")
    REORDER(k, j, i)
    // VhhhA[j + i*No + L*NoNo] * TBChh[L + k*No]; H1
    DGEMM_HOLES(VhhhC, TABhh, "N")
    REORDER(j, i, k)
    // VhhhA[k + i*No + L*NoNo] * TBChh[j + L*No]; H4
    DGEMM_HOLES(VhhhC, TABhh, "T")
    REORDER(k, i, j)
  }
  chrono["doubles:holes"].stop();

  chrono["doubles:particles"].start();
  { // Particle part =========================================================
    // TAphh[E + i*Nv + j*NoNv] * VBCph[E + k*Nv]; P0
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(i, j, k)
    // TAphh[E + i*Nv + k*NoNv] * VCBph[E + j*Nv]; P3
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(i, k, j)
    // TCphh[E + k*Nv + i*NoNv] * VABph[E + j*Nv]; P5
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(k, i, j)
    // TCphh[E + k*Nv + j*NoNv] * VBAph[E + i*Nv]; P2
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(k, j, i)
    // TBphh[E + j*Nv + i*NoNv] * VACph[E + k*Nv]; P1
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(j, i, k)
    // TBphh[E + j*Nv + k*NoNv] * VCAph[E + i*Nv]; P4
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(j, k, i)
  }
  chrono["doubles:particles"].stop();
}

int main(int argc, char ** argv){
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Timings chrono;
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

  chrono["total"].start();
  for (int it = 0; it < iterations; it++) {

    chrono["doubles"].start();
    doDoubles(No, Nv,
        VhhhC.data(), TABhh.data(),
        TAphh.data(), VBCph.data(),
        Tijk.data(),
        chrono);
    chrono["doubles"].stop();

  }
  chrono["total"].stop();


  // PRINT TIMINGS
  for (auto const& pair: chrono)
    LOG << std::setprecision(6) << std::setw(6)
        << pair.second.count() << "  "
        << pair.first
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

