CXX = mpic++
OPTIONS = -fopenmp

DEFINES += -DHAS_GCC -DHAS_LAPACK
LIBS_STATIC = -lblas -lgfortran -lquadmath
