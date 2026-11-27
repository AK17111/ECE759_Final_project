Makefile
# Build:  make        -> creates ./ldpc
# Test:   make test
# No -ffast-math: it reorders float math and would break exact matching with OpenMP/CUDA.
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall
OMPFLAGS ?= -fopenmp
# Mac (Apple clang has no OpenMP):  make OMPFLAGS=     -> builds without OpenMP (runs serial)
 
ldpc: main.cpp ldpc_code.cpp channel.cpp decoder.cpp ldpc.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) main.cpp ldpc_code.cpp channel.cpp decoder.cpp -o ldpc
 
test: ldpc
	./ldpc test
 
clean:
	rm -f ldpc
 
.PHONY: test clean