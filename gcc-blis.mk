CXX = mpic++
OPTIONS = -fopenmp
BLIS_ARCH ?= intel64
#BLIS_ARCH ?= haswell

BLIS_LIB_PATH = ./blis/lib/$(BLIS_ARCH)
LIBS_STATIC = -L$(BLIS_LIB_PATH) -lblis

DEFINES += -DHAS_GCC -DHAS_BLIS -DBLIS_ARCH=$(BLIS_ARCH)
INCLUDES += -I./blis/include/$(BLIS_ARCH)/

blis: $(BLIS_LIB_PATH)

./blis/lib/%:
	git clone https://github.com/flame/blis
	cd blis && ./configure $* && $(MAKE)
