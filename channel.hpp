// channel.hpp - test-data generation and the shared frame file format.
//
// A "frame" is one codeword. We create random payload bits, encode them,
// send them through a simulated noisy channel (BPSK over AWGN), and turn
// each received value y into a soft confidence:  LLR = 2*y / sigma^2.
//   LLR > 0  -> bit is probably 0     |LLR| large -> very confident
//   LLR < 0  -> bit is probably 1     |LLR| small -> unsure
//
// Frames are written to a binary file so the serial, OpenMP and CUDA
// decoders all decode EXACTLY the same inputs.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "ldpc_code.hpp"

struct FrameSet {
    int n = 0;                 // bits per frame
    int frames = 0;
    float ebno_db = 0.0f;      // noise level used (Eb/N0 in dB)
    std::vector<float> llr;    // frames * n channel LLRs
    std::vector<uint8_t> truth;  // frames * n transmitted bits (for error counting)
};

// sigma for a given Eb/N0 (dB) and code rate, BPSK.
double noise_sigma(double ebno_db, double rate);

FrameSet make_frames(const LdpcCode& code, int frames, double ebno_db, uint64_t seed);

// File layout (little-endian):
//   char magic[4] = "LDPF"; int32 n; int32 frames; float ebno_db;
//   float llr[frames*n]; uint8 truth[frames*n]
void write_frames(const FrameSet& fs, const std::string& path);
FrameSet read_frames(const std::string& path);

// Decoded output file: char magic[4]="LDPD"; int32 n; int32 frames;
//   uint8 bits[frames*n]; int32 iters[frames]
void write_decoded(const std::string& path, int n, int frames,
                   const std::vector<uint8_t>& bits, const std::vector<int>& iters);
void read_decoded(const std::string& path, int& n, int& frames,
                  std::vector<uint8_t>& bits, std::vector<int>& iters);
