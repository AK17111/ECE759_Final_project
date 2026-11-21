// ldpc_code.hpp - parity-check matrix H stored as a sparse Tanner graph.
//
// Every "1" in H is one EDGE of the Tanner graph. Edges are numbered in
// row (check) order, so check c owns edges chk_ptr[c] .. chk_ptr[c+1]-1.
// All decoders (serial, OpenMP, CUDA) store one message per edge, indexed
// by this same edge id, so their memory layouts line up.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct LdpcCode {
    int n = 0;   // number of bits (variable nodes, columns of H)
    int m = 0;   // number of parity checks (check nodes, rows of H)

    // Check side (CSR by rows of H)
    std::vector<int> chk_ptr;   // size m+1
    std::vector<int> edge_var;  // size E: which bit edge e connects to

    // Bit side (for each bit, the ids of its edges)
    std::vector<int> var_ptr;   // size n+1
    std::vector<int> var_edge;  // size E: edge ids, ascending per bit

    // Systematic encoder data (filled by build_encoder)
    int k = 0;                          // information bits = n - rank(H)
    std::vector<int> info_pos;          // which bit positions carry payload
    std::vector<int> parity_pos;        // pivot column for each encoder row
    std::vector<std::vector<uint64_t>> enc_rows;  // dense rows over n bits

    int num_edges() const { return (int)edge_var.size(); }
    double rate() const { return n ? (double)k / n : 0.0; }
};

// Build from a list of rows: rows[c] = bit indices in check c.
LdpcCode code_from_rows(int n, const std::vector<std::vector<int>>& rows);

// MacKay "alist" format (used by 96.33.964 and most published codes).
LdpcCode load_alist(const std::string& path);
void save_alist(const LdpcCode& code, const std::string& path);

// Progressive Edge Growth: (dv, dc)-regular code with large girth.
LdpcCode make_peg_code(int n, int dv, int dc, uint64_t seed);

// Number of length-4 cycles (pairs of bits that share 2+ checks).
long count_4cycles(const LdpcCode& code);

// Gaussian elimination over GF(2) -> systematic encoder.
void build_encoder(LdpcCode& code);
// codeword[info_pos[i]] = info[i]; parity bits filled in.
void encode(const LdpcCode& code, const uint8_t* info, uint8_t* codeword);
// true if H * c^T == 0 (mod 2)
bool syndrome_ok(const LdpcCode& code, const uint8_t* c);
