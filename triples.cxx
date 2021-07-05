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

//#undef REORDER
//#define REORDER(__II, __JJ, __KK)

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
  std::vector<double> _t_buffer;
  _t_buffer.reserve(NoNoNo);

  for (size_t k = 0; k < NoNoNo; k++) {
    // zero the Tijk
    Tijk[k] = 0.0;
  }

  chrono["doubles:holes"].start();
  { // Holes part ============================================================
    // VhhhC[i + k*No + L*NoNo] * TABhh[L + j*No]; H1
    chrono["doubles:holes:1"].start();
    DGEMM_HOLES(VhhhC, TABhh, "N")
    REORDER(i, k, j)
    chrono["doubles:holes:1"].stop();
    // VhhhC[j + k*No + L*NoNo] * TABhh[i + L*No]; H0
    chrono["doubles:holes:2"].start();
    DGEMM_HOLES(VhhhC, TABhh, "T")
    REORDER(j, k, i)
    chrono["doubles:holes:2"].stop();
    // VhhhB[i + j*No + L*NoNo] * TAChh[L + k*No]; H5
    chrono["doubles:holes:3"].start();
    DGEMM_HOLES(VhhhB, TAChh, "N")
    REORDER(i, j, k)
    chrono["doubles:holes:3"].stop();
    // VhhhB[k + j*No + L*NoNo] * TAChh[i + L*No]; H3
    chrono["doubles:holes:4"].start();
    DGEMM_HOLES(VhhhB, TAChh, "T")
    REORDER(k, j, i)
    chrono["doubles:holes:4"].stop();
    // VhhhA[j + i*No + L*NoNo] * TBChh[L + k*No]; H1
    chrono["doubles:holes:5"].start();
    DGEMM_HOLES(VhhhA, TBChh, "N")
    REORDER(j, i, k)
    chrono["doubles:holes:5"].stop();
    // VhhhA[k + i*No + L*NoNo] * TBChh[j + L*No]; H4
    chrono["doubles:holes:6"].start();
    DGEMM_HOLES(VhhhA, TBChh, "T")
    REORDER(k, i, j)
    chrono["doubles:holes:6"].stop();
  }
  chrono["doubles:holes"].stop();

  chrono["doubles:particles"].start();
  { // Particle part =========================================================
    // TAphh[E + i*Nv + j*NoNv] * VBCph[E + k*Nv]; P0
    chrono["doubles:particles:1"].start();
    DGEMM_PARTICLES(TAphh, VBCph)
    REORDER(i, j, k)
    chrono["doubles:particles:1"].stop();
    // TAphh[E + i*Nv + k*NoNv] * VCBph[E + j*Nv]; P3
    chrono["doubles:particles:2"].start();
    DGEMM_PARTICLES(TAphh, VCBph)
    REORDER(i, k, j)
    chrono["doubles:particles:2"].stop();
    // TCphh[E + k*Nv + i*NoNv] * VABph[E + j*Nv]; P5
    chrono["doubles:particles:3"].start();
    DGEMM_PARTICLES(TCphh, VABph)
    REORDER(k, i, j)
    chrono["doubles:particles:3"].stop();
    // TCphh[E + k*Nv + j*NoNv] * VBAph[E + i*Nv]; P2
    chrono["doubles:particles:4"].start();
    DGEMM_PARTICLES(TCphh, VBAph)
    REORDER(k, j, i)
    chrono["doubles:particles:4"].stop();
    // TBphh[E + j*Nv + i*NoNv] * VACph[E + k*Nv]; P1
    chrono["doubles:particles:5"].start();
    DGEMM_PARTICLES(TBphh, VACph)
    REORDER(j, i, k)
    chrono["doubles:particles:5"].stop();
    // TBphh[E + j*Nv + k*NoNv] * VCAph[E + i*Nv]; P4
    chrono["doubles:particles:6"].start();
    DGEMM_PARTICLES(TBphh, VCAph)
    REORDER(j, k, i)
    chrono["doubles:particles:6"].stop();

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

  const bool doRest(hauta::option<bool>(argc, argv, "--rest", false));

  const int NoNo = No*No, NoNoNo = No*No*No;
  const double flopCountParticles = double(NoNoNo) * double(Nv) * 2.0 * 6.0 / 1e9;
  const double flopCountHoles = double(NoNoNo) * double(No) * 2.0 * 6.0 / 1e9;
  const double flopCount = flopCountHoles + flopCountParticles;


  std::vector<double> Tijk
                    , VhhhC, TABhh
                    , VhhhB, TAChh
                    , VhhhA, TBChh
                    , TAphh, VBCph, VCBph
                    , TCphh, VABph, VBAph
                    , TBphh, VACph, VCAph
                    ;
  Tijk.reserve(NoNoNo);
  VhhhC.reserve(NoNoNo); TABhh.reserve(NoNo);
  VhhhB.reserve(NoNoNo); TAChh.reserve(NoNo);
  VhhhA.reserve(NoNoNo); TBChh.reserve(NoNo);
  TAphh.reserve(NoNo*Nv); VBCph.reserve(Nv*No); VCBph.reserve(Nv*No);
  TCphh.reserve(NoNo*Nv); VABph.reserve(Nv*No); VBAph.reserve(Nv*No);
  TBphh.reserve(NoNo*Nv); VACph.reserve(Nv*No); VCAph.reserve(Nv*No);

  std::vector<double> _Tijk
                    , _VhhhC, _TABhh
                    , _VhhhB, _TAChh
                    , _VhhhA, _TBChh
                    , _TAphh, _VBCph, _VCBph
                    , _TCphh, _VABph, _VBAph
                    , _TBphh, _VACph, _VCAph
                    ;
  _Tijk.reserve(NoNoNo);
  _VhhhC.reserve(NoNoNo);  _TABhh.reserve(NoNo);
  _VhhhB.reserve(NoNoNo);  _TAChh.reserve(NoNo);
  _VhhhA.reserve(NoNoNo);  _TBChh.reserve(NoNo);
  _TAphh.reserve(NoNo*Nv); _VBCph.reserve(Nv*No); _VCBph.reserve(Nv*No);
  _TCphh.reserve(NoNo*Nv); _VABph.reserve(Nv*No); _VBAph.reserve(Nv*No);
  _TBphh.reserve(NoNo*Nv); _VACph.reserve(Nv*No); _VCAph.reserve(Nv*No);

  LOG << "======= TRIPLES ======\n";
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

  chrono["total"].start();
  for (size_t it = 0; it < iterations; it++) {
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

    if (doRest) {
      // do something else
      chrono["rest"].start();
      doDoubles(No, Nv,
                _VhhhC.data(), _TABhh.data(),
                _VhhhB.data(), _TAChh.data(),
                _VhhhA.data(), _TBChh.data(),
                _TAphh.data(), _VBCph.data(), _VCBph.data(),
                _TCphh.data(), _VABph.data(), _VBAph.data(),
                _TBphh.data(), _VACph.data(), _VCAph.data(),
                _Tijk.data(),
                chrono);
      chrono["rest"].stop();
    }

  }
  chrono["total"].stop();


  // PRINT TIMINGS
  for (auto const& pair: chrono)
    LOG << std::setprecision(6) << std::setw(6)
        << pair.first << ": "
        << pair.second.count()
        << "\n";

  LOG
    << "flops:doubles:no-reorder "
    << flopCount * iterations
    / (chrono["doubles"].count() - chrono["doubles:reorder"].count())
    << "\n";
  LOG
    << "flops:holes:no-reorder "
    << flopCountHoles * iterations
    / ( chrono["doubles:holes"].count()
      - 0.5 * chrono["doubles:reorder"].count()
      )
   << "\n";
  LOG
    << "flops:doubles:no-reorder "
    << flopCountParticles * iterations
    / ( chrono["doubles:particles"].count()
      - 0.5 * chrono["doubles:reorder"].count()
      )
   << "\n";
  LOG
    << "flops:holes "
    << flopCountHoles * iterations / chrono["doubles:holes"].count()
    << "\n";
  LOG
    << "flops:particles "
    << flopCountParticles * iterations / chrono["doubles:particles"].count()
    << "\n";
  LOG
    << "flops:doubles "
    << flopCount * iterations / chrono["doubles"].count()
    << "\n";

  MPI_Finalize();

  return 0;

}

