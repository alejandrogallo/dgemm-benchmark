MKL = -mkl -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

CC = mpiicc
CXX = mpiicc
OPTIONS = \
-fopenmp $(MKL) -qoverride-limits \
-Qoption,cpp,--extended_float_types


DEFINES += -DHAS_INTEL -DHAS_MKL
