MKL = -mkl -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

CXX = mpic++
OPTIONS = \
-fopenmp $(MKL) -qoverride-limits \
-qopenmp -Qoption,cpp,--extended_float_types


DEFINES += -DHAS_INTEL -DHAS_MKL
