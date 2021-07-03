SRC = dgemm.cxx
TARGETS = dgemm-debug dgemm

CONFIG ?= gcc
include ${CONFIG}.mk


all: $(TARGETS)
$(TARGETS): $(SRC) Makefile
CXX_FLAGS =   \
-std=c++11    \
$(OPTIONS)    \
-march=native \
-Wl,-Bstatic $(LIBS_STATIC) -Wl,-Bdynamic


%-debug: %.cxx
	$(CXX) $(DEFINES) $(INCLUDES) -DDEBUG='"yes"' $< $(CXX_FLAGS) -O0 -g -o $@

%: %.cxx
	$(CXX) $(DEFINES) $(INCLUDES) -DDEBUG='"no"' $< $(CXX_FLAGS) -O3 -o $@

clean:
	rm -f $(TARGETS)

.PHONY: all clean
