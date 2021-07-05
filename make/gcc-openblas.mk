CXX = mpic++
OPTIONS = -fopenmp

DEFINES += -DHAS_GCC -DHAS_OPENBLAS
LIBS_STATIC = -lopenblas -lquadmath
