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

#undef REORDER
#define REORDER(__II, __JJ, __KK)

void doDoubles(size_t No, size_t Nv,
               double const* VhhhC, double const* TABhh,
               double const* VhhhB, double const* TAChh,
               double const* VhhhA, double const* TBChh,
               // particles
               double const* TAphh, double const* VBCph, double const* VCBph,
               double const* TCphh, double const* VABph, double const* VBAph,
               double const* TBphh, double const* VACph, double const* VCAph,
               double* Tijk,
               Timings& chrono) {

  const size_t NoNo = No*No, NoNoNo = No*No*No;
  double one(1.0), m_one(-1.0), zero(0.0);
  std::vector<double> _t_buffer(NoNoNo);

  for (size_t k = 0; k < NoNoNo; k++) {
    // zero the Tijk
    Tijk[k] = 0.0;
  }

  chrono["doubles:holes"].start();
  { // Holes part ============================================================
    // VhhhC[i + k*No + L*NoNo] * TABhh[L + j*No]; H1
    DGEMM_HOLES(VhhhC, TABhh, "N")
    REORDER(i, k, j)
    // VhhhC[j + k*No + L*NoNo] * TABhh[i + L*No]; H0
    DGEMM_HOLES(VhhhC, TABhh, "T")
    REORDER(j, k, i)
    // VhhhB[i + j*No + L*NoNo] * TAChh[L + k*No]; H5
    DGEMM_HOLES(VhhhB, TAChh, "N")
    REORDER(i, j, k)
    // VhhhB[k + j*No + L*NoNo] * TAChh[i + L*No]; H3
    DGEMM_HOLES(VhhhB, TAChh, "T")
    REORDER(k, j, i)
    // VhhhA[j + i*No + L*NoNo] * TBChh[L + k*No]; H1
    DGEMM_HOLES(VhhhA, TBChh, "N")
    REORDER(j, i, k)
    // VhhhA[k + i*No + L*NoNo] * TBChh[j + L*No]; H4
    DGEMM_HOLES(VhhhA, TBChh, "T")
    REORDER(k, i, j)
  }
  chrono["doubles:holes"].stop();

  chrono["doubles:particles"].start();
  { // Particle part =========================================================
    // TAphh[E + i*Nv + j*NoNv] * VBCph[E + k*Nv]; P0
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(i, j, k)
    // TAphh[E + i*Nv + k*NoNv] * VCBph[E + j*Nv]; P3
    DGEMM_PARTICLES(TAphh, VCBph)
    REORDER(i, k, j)
    // TCphh[E + k*Nv + i*NoNv] * VABph[E + j*Nv]; P5
    DGEMM_PARTICLES(TCphh, VABph)
    REORDER(k, i, j)
    // TCphh[E + k*Nv + j*NoNv] * VBAph[E + i*Nv]; P2
    DGEMM_PARTICLES(TCphh, VBAph)
    REORDER(k, j, i)
    // TBphh[E + j*Nv + i*NoNv] * VACph[E + k*Nv]; P1
    DGEMM_PARTICLES(TBphh, VACph)
    REORDER(j, i, k)
    // TBphh[E + j*Nv + k*NoNv] * VCAph[E + i*Nv]; P4
    DGEMM_PARTICLES(TBphh, VCAph)
    REORDER(j, k, i)

  }
  chrono["doubles:particles"].stop();
}

int main(int argc, char ** argv){
  int rank, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Timings chrono;

  size_t No = hauta::option<size_t>(argc, argv, "--no")
       , Nv = hauta::option<size_t>(argc, argv, "--nv")
       , iterations = hauta::option<size_t>(argc, argv, "-i", 1)
       ;

  const int NoNo = No*No, NoNoNo = No*No*No;
  const double flopCountParticles = double(NoNoNo) * double(Nv) * 2.0 * 6.0 / 1e9;
  const double flopCountHoles = double(NoNoNo) * double(No) * 2.0 * 6.0 / 1e9;
  const double flopCount = flopCountHoles + flopCountParticles;

  std::vector<double> Tijk(NoNoNo);
  std::vector<double>
      VhhhC(NoNoNo), TABhh(NoNo)
    , VhhhB(NoNoNo), TAChh(NoNo)
    , VhhhA(NoNoNo), TBChh(NoNo)
    ;
  std::vector<double>
      TAphh(NoNo*Nv), VBCph(Nv*No), VCBph(Nv*No)
    , TCphh(NoNo*Nv), VABph(Nv*No), VBAph(Nv*No)
    , TBphh(NoNo*Nv), VACph(Nv*No), VCAph(Nv*No)
    ;

  LOG << "Doing " << argv[0] << " Tests\n";
  LOG << "»»»»»»»»»»»»»»»»»\n";
  LOG << SHOW_VAR(np) << "\n";
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
    chrono["start:stop"].start();
    chrono["start:stop"].stop();

    chrono["doubles"].start();
    doDoubles(No, Nv,
        VhhhC.data(), TABhh.data(),
        VhhhB.data(), TAChh.data(),
        VhhhA.data(), TBChh.data(),
        TAphh.data(), VBCph.data(), VCBph.data(),
        TCphh.data(), VABph.data(), VBAph.data(),
        TBphh.data(), VACph.data(), VCAph.data(),
        Tijk.data(),
        chrono);
    chrono["doubles"].stop();

  }
  chrono["total"].stop();


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
    / (chrono["doubles"].count() - chrono["doubles:reorder"].count())
    << "\n";
  LOG
    << "flops:holes "
    << flopCountHoles * iterations / chrono["doubles:holes"].count()
    << "\n";
  LOG
    << "flops:holes:no-reorder "
    << flopCountHoles * iterations
    / ( chrono["doubles:holes"].count()
      - 0.5 * chrono["doubles:reorder"].count()
      )
   << "\n";
  LOG
    << "flops:particles "
    << flopCountParticles * iterations / chrono["doubles:particles"].count()
    << "\n";
  LOG
    << "flops:doubles:no-reorder "
    << flopCountParticles * iterations
    / ( chrono["doubles:particles"].count()
      - 0.5 * chrono["doubles:reorder"].count()
      )
   << "\n";

  MPI_Finalize();

  return 0;

}

