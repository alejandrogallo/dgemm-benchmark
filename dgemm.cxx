#include <iostream>
#include <vector>
#include <iomanip>
#include "utils.h"
#ifdef HAS_BLIS
#include <blis.h>
#endif

#if !defined(DEBUG)
#  define DEBUG "unknown"
#endif


int main(int argc, char ** argv){

  Timings chrono;
  double one(1.0), m_one(-1.0), zero(0.0);
  const int
      No(atoi(argv[1]))
    , Nv(atoi(argv[2]))
    , iterations(atoi(argv[3]))
    ;

  const int NoNo = No*No, NoNoNo = No*No*No;
  const double flopCount = double(NoNoNo) * double(No+Nv) * 2.0 * 6.0 / 1e9;

  std::vector<double>
      TABhh(No*No), VhhhC(NoNoNo), Tijk(NoNoNo)
    , TAphh(NoNo*Nv), VBCph(Nv*No);

  std::cout << "Doing DGEMM Tests\n";
  std::cout << "»»»»»»»»»»»»»»»»»\n";
  std::cout << "No: " << No << "\n";
  std::cout << "Nv: " << Nv << "\n";
  std::cout << "its: " << iterations << "\n";
  std::cout << "debug: " << DEBUG << "\n";
#if defined(HAS_INTEL)
  std::cout << "intel compiler\n";
#endif

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
      chrono["holes:reorder"].start();
      for (size_t ijk = 0; ijk < NoNoNo; ijk++) {
        Tijk.data()[ijk] = 1.0;
      }
      chrono["holes:reorder"].stop();
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
      chrono["particles:reorder"].start();
      for (size_t ijk = 0; ijk < NoNoNo; ijk++) {
        Tijk.data()[ijk] = 1.0;
      }
      chrono["particles:reorder"].stop();
    }
    chrono["particles"].stop();

  }
  chrono["doubles"].stop();




  // PRINT TIMINGS
  for (auto const& pair: chrono)
    std::cout << std::setprecision(6) << std::setw(6)
              << pair.second.count() << "  "
              << pair.first
              << std::endl;

  std::cout
    << "flops:doubles "
    << flopCount * iterations / chrono["doubles"].count()
    << "\n";
  std::cout
    << "flops:doubles:no-reorder "
    << flopCount * iterations
    / (chrono["particles:dgemm"].count() + chrono["holes:dgemm"].count())
    << "\n";

  return 0;

}

