# Build:  make        -> creates ./ldpc
# Test:   make test
# No -ffast-math: it reorders float math and would break exact matching with OpenMP/CUDA.
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall

ldpc: main.cpp ldpc_code.cpp channel.cpp decoder.cpp ldpc.hpp
	$(CXX) $(CXXFLAGS) main.cpp ldpc_code.cpp channel.cpp decoder.cpp -o ldpc

test: ldpc
	./ldpc test

clean:
	rm -f ldpc

.PHONY: test clean