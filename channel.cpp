// channel.cpp - noisy-frame generation and binary file I/O.
#include "channel.hpp"

#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <stdexcept>

double noise_sigma(double ebno_db, double rate) {
    double ebno = std::pow(10.0, ebno_db / 10.0);
    return std::sqrt(1.0 / (2.0 * rate * ebno));
}

FrameSet make_frames(const LdpcCode& code, int frames, double ebno_db, uint64_t seed) {
    FrameSet fs;
    fs.n = code.n;
    fs.frames = frames;
    fs.ebno_db = (float)ebno_db;
    fs.llr.resize((size_t)frames * code.n);
    fs.truth.resize((size_t)frames * code.n);
    const double sigma = noise_sigma(ebno_db, code.rate());
    const double scale = 2.0 / (sigma * sigma);
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> noise(0.0, sigma);
    std::vector<uint8_t> info(code.k);
    for (int f = 0; f < frames; ++f) {
        for (auto& b : info) b = (uint8_t)(rng() & 1);
        uint8_t* cw = &fs.truth[(size_t)f * code.n];
        encode(code, info.data(), cw);
        for (int j = 0; j < code.n; ++j) {
            double x = cw[j] ? -1.0 : +1.0;      // BPSK: bit 0 -> +1, bit 1 -> -1
            double y = x + noise(rng);           // channel adds Gaussian noise
            fs.llr[(size_t)f * code.n + j] = (float)(scale * y);
        }
    }
    return fs;
}

template <class T> static void put(std::ofstream& o, const T& v) { o.write((const char*)&v, sizeof(T)); }
template <class T> static void get(std::ifstream& i, T& v) { i.read((char*)&v, sizeof(T)); }

void write_frames(const FrameSet& fs, const std::string& path) {
    std::ofstream o(path, std::ios::binary);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.write("LDPF", 4);
    put(o, (int32_t)fs.n); put(o, (int32_t)fs.frames); put(o, fs.ebno_db);
    o.write((const char*)fs.llr.data(), fs.llr.size() * sizeof(float));
    o.write((const char*)fs.truth.data(), fs.truth.size());
}

FrameSet read_frames(const std::string& path) {
    std::ifstream i(path, std::ios::binary);
    if (!i) throw std::runtime_error("cannot open " + path);
    char magic[4];
    i.read(magic, 4);
    if (std::memcmp(magic, "LDPF", 4)) throw std::runtime_error("not a frame file: " + path);
    FrameSet fs;
    int32_t n, f;
    get(i, n); get(i, f); get(i, fs.ebno_db);
    fs.n = n; fs.frames = f;
    fs.llr.resize((size_t)n * f);
    fs.truth.resize((size_t)n * f);
    i.read((char*)fs.llr.data(), fs.llr.size() * sizeof(float));
    i.read((char*)fs.truth.data(), fs.truth.size());
    if (!i) throw std::runtime_error("truncated frame file: " + path);
    return fs;
}

void write_decoded(const std::string& path, int n, int frames,
                   const std::vector<uint8_t>& bits, const std::vector<int>& iters) {
    std::ofstream o(path, std::ios::binary);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.write("LDPD", 4);
    put(o, (int32_t)n); put(o, (int32_t)frames);
    o.write((const char*)bits.data(), bits.size());
    for (int it : iters) put(o, (int32_t)it);
}

void read_decoded(const std::string& path, int& n, int& frames,
                  std::vector<uint8_t>& bits, std::vector<int>& iters) {
    std::ifstream i(path, std::ios::binary);
    if (!i) throw std::runtime_error("cannot open " + path);
    char magic[4];
    i.read(magic, 4);
    if (std::memcmp(magic, "LDPD", 4)) throw std::runtime_error("not a decoded file: " + path);
    int32_t n32, f32;
    get(i, n32); get(i, f32);
    n = n32; frames = f32;
    bits.resize((size_t)n * frames);
    i.read((char*)bits.data(), bits.size());
    iters.resize(frames);
    for (auto& it : iters) { int32_t v; get(i, v); it = v; }
    if (!i) throw std::runtime_error("truncated decoded file: " + path);
}
