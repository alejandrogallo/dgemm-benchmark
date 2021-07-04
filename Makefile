SRC = dgemm.cxx
TARGETS = dgemm-debug dgemm
TARGETS += triples-debug triples

CONFIG ?= icc
include ${CONFIG}.mk

$(info CONFIG = $(CONFIG))


all: $(patsubst %,bin/${CONFIG}/%,$(TARGETS))
$(TARGETS): $(SRC) Makefile
CXX_FLAGS =     \
-pedantic -Wall \
-std=c++11      \
-fmax-errors=1  \
$(OPTIONS)      \
-march=native   \
-Wl,-Bstatic $(LIBS_STATIC) -Wl,-Bdynamic

DEFINES += -DGIT_COMMIT="$(shell git describe --always)"
DEFINES += -DCONFIG=$(CONFIG)

bin/${CONFIG}/%-debug: %.cxx
	@mkdir -p ${@D}
	$(CXX) $(DEFINES) $(INCLUDES) -DDEBUG $< $(CXX_FLAGS) -O0 -g -o $@

bin/${CONFIG}/%: %.cxx
	@mkdir -p ${@D}
	$(CXX) $(DEFINES) $(INCLUDES) $< $(CXX_FLAGS) -O3 -o $@

clean:
	rm -r bin

.PHONY: all clean
