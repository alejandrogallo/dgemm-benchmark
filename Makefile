-include config.mk
.DEFAULT_GOAL := bin

TARGETS += dgemm-debug   dgemm
TARGETS += triples-debug triples
TARGETS += blas-debug    blas
TARGETS += vector-debug  vector

BIN := $(patsubst %,bin/${CONFIG}/%,$(TARGETS))
OBJ := $(patsubst %,obj/${CONFIG}/%.o,$(TARGETS))
LST := $(patsubst %,obj/${CONFIG}/%.lst,$(TARGETS))

CONFIG ?= icc
include make/${CONFIG}.mk
$(info CONFIG = $(CONFIG))


CXX_FLAGS =     \
-pedantic -Wall \
-ansi \
-std=c++11      \
-fmax-errors=1  \
$(OPTIONS)      \
-march=native   \
-Wl,-Bstatic $(LIBS_STATIC) -Wl,-Bdynamic

CFLAGS =        \
-pedantic -Wall \
-ansi \
-fmax-errors=1  \
$(OPTIONS)      \
-march=native   \
-Wl,-Bstatic $(LIBS_STATIC) -Wl,-Bdynamic

DEB_FLAGS = -DDEBU -O0 -g
OPT_FLAGS = -O3

all: bin obj lst
bin: $(BIN)
obj: $(OBJ)
lst: $(LST) Makefile


$(TARGETS): Makefile

DEFINES += -DGIT_COMMIT="$(shell git describe --always)"
DEFINES += -DCONFIG=$(CONFIG)
DEFINES += -DCOMPILER_VERSION="$(shell $(CXX) --version)"
DEFINES += -DDATE="$(shell date)"

bin/${CONFIG}/%-debug: %.cxx
	@mkdir -p ${@D}
	$(CXX) $(DEFINES) $(INCLUDES) ${DEB_FLAGS} $< $(CXX_FLAGS) -o $@

bin/${CONFIG}/%: %.cxx
	@mkdir -p ${@D}
	$(CXX) $(DEFINES) $(INCLUDES) ${OPT_FLAGS} $< $(CXX_FLAGS) -o $@

bin/${CONFIG}/%-debug: %.c
	@mkdir -p ${@D}
	$(CC) $(DEFINES) $(INCLUDES) ${DEB_FLAGS} $< $(CFLAGS) -o $@

bin/${CONFIG}/%: %.c
	@mkdir -p ${@D}
	$(CC) $(DEFINES) $(INCLUDES) ${OPT_FLAGS} $< $(CFLAGS) -o $@


obj/${CONFIG}/%-debug.o: %.cxx
	@mkdir -p ${@D}
	$(CXX) $(DEFINES) $(INCLUDES) ${DEB_FLAGS} -c $< $(CXX_FLAGS) -o $@

obj/${CONFIG}/%.o: %.cxx
	@mkdir -p ${@D}
	$(CXX) $(DEFINES) $(INCLUDES) ${OPT_FLAGS} -c $< $(CXX_FLAGS) -o $@

%.lst: %.o
	objdump --demangle -f -l -Mintel -S $< > $@

clean:
	rm -r bin/${CONFIG} obj/${CONFIG}

.PHONY: all clean bin obj lst
