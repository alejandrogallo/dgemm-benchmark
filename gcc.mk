CXX = mpic++
OPTIONS = -fopenmp

LIBS_STATIC = -Lblis/lib/haswell/ -lblis

DEFINES += -DHAS_GCC -DHAS_BLIS
INCLUDES += -I./blis/include/haswell/

blis: ./blis/lib/haswell/libblis.a

./blis/lib/%/libblis.a:
	git clone https://github.com/flame/blis
	cd blis && ./configure $* && $(MAKE)
