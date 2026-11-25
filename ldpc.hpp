#pragma once
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

//  H matrix
struct LdpcCode {
    int n = 0;   // number of bits    (columns of H)
    int m = 0;   // number of checks  (rows of H)

    std::vector<int> chk_ptr;   // size m+1   check side
    std::vector<int> edge_var;  // size E     edge e connects to bit edge_var[e]
    std::vector<int> var_ptr;   // size n+1   bit side
    std::vector<int> var_edge;  // size E     bit j's edge ids: var_edge[var_ptr[j] .. var_ptr[j+1]-1]

    // Encoder data (filled by build_encoder)
    int k = 0;                                    // payload bits
    std::vector<int> info_pos;                    // where payload bits sit in the codeword
    std::vector<int> parity_pos;                  // where each parity bit sits
    std::vector<std::vector<uint64_t>> enc_rows;  // which payload bits each parity bit XORs

    int num_edges() const { return (int)edge_var.size(); }
    double rate() const { return n ? (double)k / n : 0.0; }
};

LdpcCode code_from_rows(int n, const std::vector<std::vector<int>>& rows);
LdpcCode load_alist(const std::string& path);
void save_alist(const LdpcCode& code, const std::string& path);
LdpcCode make_peg_code(int n, int dv, int dc, uint64_t seed);
long count_4cycles(const LdpcCode& code);
void build_encoder(LdpcCode& code);
void encode(const LdpcCode& code, const uint8_t* info, uint8_t* codeword);
bool syndrome_ok(const LdpcCode& code, const uint8_t* c);

//noise + files
struct FrameSet {
    int n = 0;                   // bits per codeword
    int frames = 0;              // number of codewords
    float ebno_db = 0.0f;        // noise level used (Eb/N0 in dB)
    std::vector<float> llr;      // frames*n  decoder input  (LLR = 2y/sigma^2)
    std::vector<uint8_t> truth;  // frames*n  original bits  (to count errors)
};

double noise_sigma(double ebno_db, double rate);
FrameSet make_frames(const LdpcCode& code, int frames, double ebno_db, uint64_t seed);
void write_frames(const FrameSet& fs, const std::string& path);
FrameSet read_frames(const std::string& path);
void write_decoded(const std::string& path, int n, int frames,
                   const std::vector<uint8_t>& bits, const std::vector<int>& iters);
void read_decoded(const std::string& path, int& n, int& frames,
                  std::vector<uint8_t>& bits, std::vector<int>& iters);

//  decoder
struct DecodeParams {
    int max_iters = 20;
    float alpha = 0.75f;      // scaling factor (1.0 = plain min-sum)
    bool early_stop = true;   // stop as soon as every parity check is satisfied
};

// Scratch memory: one message per edge in each direction. Reused for every codeword.
// (OpenMP: every thread needs its OWN workspace.)
struct DecoderWorkspace {
    std::vector<float> v2c;   // bit   -> check messages
    std::vector<float> c2v;   // check -> bit   messages
    explicit DecoderWorkspace(const LdpcCode& code)
        : v2c(code.num_edges()), c2v(code.num_edges()) {}
};

int decode_frame(const LdpcCode& code, const float* llr, uint8_t* bits,
                 const DecodeParams& p, DecoderWorkspace& ws, bool* ok);
void decode_batch_serial(const LdpcCode& code, int frames, const float* llr,
                         uint8_t* bits, int* iters, const DecodeParams& p);

//  command-line helper
struct Args {
    int argc; char** argv;
    Args(int c, char** v) : argc(c), argv(v) {}
    const char* find(const char* name) const {
        for (int i = 1; i < argc; ++i)
            if (!std::strcmp(argv[i], name)) return i + 1 < argc ? argv[i + 1] : "";
        return nullptr;
    }
    bool has(const char* name) const { return find(name) != nullptr; }
    std::string str(const char* name, const std::string& def = "") const {
        const char* s = find(name); return s ? s : def;
    }
    long num(const char* name, long def) const { const char* s = find(name); return s ? std::atol(s) : def; }
    double real(const char* name, double def) const { const char* s = find(name); return s ? std::atof(s) : def; }
    std::string need(const char* name) const {
        const char* s = find(name);
        if (!s || !*s) throw std::runtime_error(std::string("missing required option ") + name);
        return s;
    }
};