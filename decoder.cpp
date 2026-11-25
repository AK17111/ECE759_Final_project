// decoder.cpp  flooding normalized min-sum

//   Step 1  init:        v2c[e] = Lch(bit of e)
//   Step 2  check node:  c2v[e] = alpha * (sign of the OTHER inputs) * (min |OTHER inputs|)
//   Step 3  bit node:    v2c[e] = Lch + sum of the OTHER incoming c2v
//   Step 4  decision:    total = Lch + sum of ALL incoming c2v; bit = (total < 0)
//   Step 5  syndrome:    stop if every check is satisfied

#include "ldpc.hpp"

#include <cmath>
#include <limits>

static bool all_checks_satisfied(const LdpcCode& code, const uint8_t* bits) {
    for (int c = 0; c < code.m; ++c) {
        uint8_t parity = 0;
        for (int e = code.chk_ptr[c]; e < code.chk_ptr[c + 1]; ++e) parity ^= bits[code.edge_var[e]];
        if (parity) return false;
    }
    return true;
}

int decode_frame(const LdpcCode& code, const float* llr, uint8_t* bits,
                 const DecodeParams& p, DecoderWorkspace& ws, bool* ok) {
    const int E = code.num_edges();
    float* v2c = ws.v2c.data();
    float* c2v = ws.c2v.data();

    // Step 1: every edge starts with its bit's channel LLR.
    for (int e = 0; e < E; ++e) v2c[e] = llr[code.edge_var[e]];
    for (int j = 0; j < code.n; ++j) bits[j] = llr[j] < 0.0f;
    if (p.early_stop && all_checks_satisfied(code, bits)) { *ok = true; return 0; }

    int it = 0;
    for (it = 1; it <= p.max_iters; ++it) {
        // Step 2: check-node update. Only the two smallest magnitudes are needed:
        // the edge that holds the smallest gets the second smallest, everyone else
        // gets the smallest. The sign is the product of the others' signs.
        for (int c = 0; c < code.m; ++c) {
            const int beg = code.chk_ptr[c], end = code.chk_ptr[c + 1];
            float min1 = std::numeric_limits<float>::max(), min2 = min1;
            int min_edge = -1;
            int neg = 0;  // parity of the number of negative inputs
            for (int e = beg; e < end; ++e) {
                const float x = v2c[e];
                const float a = std::fabs(x);
                neg ^= (x < 0.0f);
                if (a < min1) { min2 = min1; min1 = a; min_edge = e; }
                else if (a < min2) { min2 = a; }
            }
            for (int e = beg; e < end; ++e) {
                const float mag = p.alpha * (e == min_edge ? min2 : min1);
                const int s = neg ^ (v2c[e] < 0.0f);  // remove this edge's own sign
                c2v[e] = s ? -mag : mag;
            }
        }
        // Steps 3 + 4: bit-node update and hard decision.
        for (int j = 0; j < code.n; ++j) {
            const int beg = code.var_ptr[j], end = code.var_ptr[j + 1];
            float total = llr[j];
            for (int k = beg; k < end; ++k) total += c2v[code.var_edge[k]];
            bits[j] = total < 0.0f;
            for (int k = beg; k < end; ++k) {
                const int e = code.var_edge[k];
                v2c[e] = total - c2v[e];  // "all" minus "own" = extrinsic
            }
        }
        // Step 5: syndrome check / early stopping.
        if (p.early_stop && all_checks_satisfied(code, bits)) { *ok = true; return it; }
    }
    *ok = all_checks_satisfied(code, bits);
    return p.max_iters;
}

void decode_batch_serial(const LdpcCode& code, int frames, const float* llr,
                         uint8_t* bits, int* iters, const DecodeParams& p) {
    DecoderWorkspace ws(code);
    for (int f = 0; f < frames; ++f) {
        bool ok;
        iters[f] = decode_frame(code, llr + (size_t)f * code.n, bits + (size_t)f * code.n, p, ws, &ok);
    }
}