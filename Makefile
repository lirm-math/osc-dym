CXX      = g++
CXXFLAGS = -fmax-errors=1 -std=c++20 -march=x86-64 -mtune=generic -I.
LDFLAGS  = -L.

REAL_4_FLAGS  = -DREAL_4
REAL_4_LIBS   = -lsolve_dym_4

REAL_16_FLAGS = -DREAL_16
REAL_16_LIBS  = -lsolve_dym_16

REAL_10_FLAGS = -DREAL_10
REAL_10_LIBS  = -lsolve_dym_10

REAL_8_FLAGS  = -DREAL_8
REAL_8_LIBS   = -lsolve_dym_8

BASIC_TARGETS = out/osc_dym_4 out/osc_dym_16 out/osc_dym_10 out/osc_dym_8

MORE_TARGETS  = out/osc_dym_more_4 out/osc_dym_more_16 \
                out/osc_dym_more_10 out/osc_dym_more_8

basic: $(BASIC_TARGETS)

more: $(MORE_TARGETS)

all: basic more

out/%_4: %.cpp
	mkdir -p out
	$(CXX) $(CXXFLAGS) $(REAL_4_FLAGS) $< $(LDFLAGS) $(REAL_4_LIBS) -o $@

out/%_16: %.cpp
	mkdir -p out
	$(CXX) $(CXXFLAGS) $(REAL_16_FLAGS) $< $(LDFLAGS) $(REAL_16_LIBS) -o $@

out/%_10: %.cpp
	mkdir -p out
	$(CXX) $(CXXFLAGS) $(REAL_10_FLAGS) $< $(LDFLAGS) $(REAL_10_LIBS) -o $@

out/%_8: %.cpp
	mkdir -p out
	$(CXX) $(CXXFLAGS) $(REAL_8_FLAGS) $< $(LDFLAGS) $(REAL_8_LIBS) -o $@

clean:
	rm -rf out

.PHONY: basic more all clean