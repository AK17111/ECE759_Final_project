// ldpc_code.cpp- THE H MATRIX and THE ENCODER: build/load/save codes, encode payloads.
#include "ldpc.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>

LdpcCode code_from_rows(int n, const std::vector<std::vector<int>>& rows) {
    LdpcCode code;
    code.n = n;                 // number of bits  (columns of H)
    code.m = (int)rows.size();  // number of rules (rows of H)

    code.chk_ptr.push_back(0);
    for (int c = 0; c < code.m; c++) {
        std::vector<int> bits = rows[c];
        std::sort(bits.begin(), bits.end());                 // keep each rule's bits in order
        for (int i = 0; i < (int)bits.size(); i++) {
            if (i > 0 && bits[i] == bits[i - 1]) continue;   // skip a bit listed twice
            if (bits[i] < 0 || bits[i] >= n) throw std::runtime_error("bit index out of range");
            code.edge_var.push_back(bits[i]);                // add one edge
        }
        code.chk_ptr.push_back((int)code.edge_var.size()); // next rule starts here
    }


    std::vector<std::vector<int>> edges_of_bit(n);
    for (int e = 0; e < code.num_edges(); e++) {
        int bit = code.edge_var[e];
        edges_of_bit[bit].push_back(e);
    }

    code.var_ptr.push_back(0);
    for (int bit = 0; bit < n; bit++) {
        for (int e : edges_of_bit[bit]) code.var_edge.push_back(e);
        code.var_ptr.push_back((int)code.var_edge.size());
    }
    return code;
}

//  alist I/O
LdpcCode load_alist(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    int n, m, max_col, max_row;
    in >> n >> m >> max_col >> max_row;
    std::vector<int> col_w(n), row_w(m);
    for (int& w : col_w) in >> w;
    for (int& w : row_w) in >> w;
    // Column section: 1-based check indices; 0 is padding and is skipped.
    std::vector<std::vector<int>> rows(m);
    for (int v = 0; v < n; ++v) {
        int got = 0;
        while (got < col_w[v]) {
            int c;
            if (!(in >> c)) throw std::runtime_error("truncated alist file");
            if (c == 0) continue;
            rows[c - 1].push_back(v);
            ++got;
        }
    }
    return code_from_rows(n, rows);
}

void save_alist(const LdpcCode& code, const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path + " (does the folder exist?)");
    std::vector<std::vector<int>> cols(code.n);
    int max_row = 0, max_col = 0;
    for (int c = 0; c < code.m; ++c) {
        max_row = std::max(max_row, code.chk_ptr[c + 1] - code.chk_ptr[c]);
        for (int e = code.chk_ptr[c]; e < code.chk_ptr[c + 1]; ++e)
            cols[code.edge_var[e]].push_back(c);
    }
    for (auto& col : cols) max_col = std::max(max_col, (int)col.size());
    out << code.n << " " << code.m << "\n" << max_col << " " << max_row << "\n";
    for (auto& col : cols) out << col.size() << " ";
    out << "\n";
    for (int c = 0; c < code.m; ++c) out << code.chk_ptr[c + 1] - code.chk_ptr[c] << " ";
    out << "\n";
    for (auto& col : cols) {
        for (int i = 0; i < max_col; ++i) out << (i < (int)col.size() ? col[i] + 1 : 0) << " ";
        out << "\n";
    }
    for (int c = 0; c < code.m; ++c) {
        int deg = code.chk_ptr[c + 1] - code.chk_ptr[c];
        for (int i = 0; i < max_row; ++i)
            out << (i < deg ? code.edge_var[code.chk_ptr[c] + i] + 1 : 0) << " ";
        out << "\n";
    }
}

// PEG construction
// For each bit, add its dv edges one at a time. Each new edge goes to the
// check that is FARTHEST from the bit in the graph built so far (so no short
// cycles are created), breaking ties by lowest current check degree.
LdpcCode make_peg_code(int n, int dv, int dc, uint64_t seed) {
    if ((long)n * dv % dc != 0) throw std::runtime_error("n*dv must be divisible by dc");
    const int m = (int)((long)n * dv / dc);
    std::mt19937_64 rng(seed);
    std::vector<std::vector<int>> chk_nbrs(m), var_nbrs(n);
    std::vector<int> dist_c(m), dist_v(n);
    const int INF = std::numeric_limits<int>::max();

    for (int v = 0; v < n; ++v) {
        for (int k = 0; k < dv; ++k) {
            // BFS from v over the current graph; dist_c = depth at which check is reached.
            std::fill(dist_c.begin(), dist_c.end(), INF);
            std::fill(dist_v.begin(), dist_v.end(), INF);
            std::queue<int> q;  // holds variable nodes
            dist_v[v] = 0;
            q.push(v);
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int c : var_nbrs[u]) {
                    if (dist_c[c] != INF) continue;
                    dist_c[c] = dist_v[u];
                    for (int w : chk_nbrs[c])
                        if (dist_v[w] == INF) { dist_v[w] = dist_v[u] + 1; q.push(w); }
                }
            }
            // Candidates: checks with room left that v is not already on.
            int best_dist = -1, best_deg = INF;
            std::vector<int> cand;
            for (int c = 0; c < m; ++c) {
                if ((int)chk_nbrs[c].size() >= dc) continue;
                if (dist_c[c] == 0) continue;  // already connected to v
                int d = dist_c[c], deg = (int)chk_nbrs[c].size();
                if (d > best_dist || (d == best_dist && deg < best_deg)) {
                    best_dist = d; best_deg = deg; cand.clear();
                }
                if (d == best_dist && deg == best_deg) cand.push_back(c);
            }
            if (cand.empty()) throw std::runtime_error("PEG ran out of check sockets");
            int c = cand[rng() % cand.size()];
            chk_nbrs[c].push_back(v);
            var_nbrs[v].push_back(c);
        }
    }
    return code_from_rows(n, chk_nbrs);
}

long count_4cycles(const LdpcCode& code) {
    // A 4-cycle = the same pair of bits appearing together in two checks.
    std::set<std::pair<int, int>> seen;
    long cycles = 0;
    for (int c = 0; c < code.m; ++c)
        for (int a = code.chk_ptr[c]; a < code.chk_ptr[c + 1]; ++a)
            for (int b = a + 1; b < code.chk_ptr[c + 1]; ++b)
                if (!seen.insert({code.edge_var[a], code.edge_var[b]}).second) ++cycles;
    return cycles;
}

// systematic encoder
void build_encoder(LdpcCode& code) {
    const int n = code.n, words = (n + 63) / 64;
    std::vector<std::vector<uint64_t>> h(code.m, std::vector<uint64_t>(words, 0));
    for (int c = 0; c < code.m; ++c)
        for (int e = code.chk_ptr[c]; e < code.chk_ptr[c + 1]; ++e) {
            int v = code.edge_var[e];
            h[c][v / 64] |= 1ULL << (v % 64);
        }
    // Reduced row echelon form; pivots searched from the last column down,
    // so parity bits land at the end and payload bits are mostly first.
    std::vector<int> pivot_col;
    std::vector<char> is_pivot(n, 0);
    int r = 0;
    for (int col = n - 1; col >= 0 && r < code.m; --col) {
        const int w = col / 64;
        const uint64_t bit = 1ULL << (col % 64);
        int sel = -1;
        for (int i = r; i < code.m; ++i) if (h[i][w] & bit) { sel = i; break; }
        if (sel < 0) continue;
        std::swap(h[r], h[sel]);
        for (int i = 0; i < code.m; ++i)
            if (i != r && (h[i][w] & bit))
                for (int x = 0; x < words; ++x) h[i][x] ^= h[r][x];
        pivot_col.push_back(col);
        is_pivot[col] = 1;
        ++r;
    }
    code.k = n - r;
    code.info_pos.clear();
    for (int v = 0; v < n; ++v) if (!is_pivot[v]) code.info_pos.push_back(v);
    code.parity_pos = pivot_col;
    code.enc_rows.assign(r, {});
    for (int i = 0; i < r; ++i) {
        h[i][pivot_col[i] / 64] &= ~(1ULL << (pivot_col[i] % 64));  // leave only info bits
        code.enc_rows[i] = std::move(h[i]);
    }
}

void encode(const LdpcCode& code, const uint8_t* info, uint8_t* cw) {
    const int words = (code.n + 63) / 64;
    std::vector<uint64_t> bits(words, 0);
    std::fill(cw, cw + code.n, 0);
    for (int i = 0; i < code.k; ++i) {
        int v = code.info_pos[i];
        cw[v] = info[i] & 1;
        if (cw[v]) bits[v / 64] |= 1ULL << (v % 64);
    }
    for (size_t i = 0; i < code.enc_rows.size(); ++i) {
        int p = 0;
        for (int x = 0; x < words; ++x) p ^= __builtin_popcountll(code.enc_rows[i][x] & bits[x]) & 1;
        cw[code.parity_pos[i]] = (uint8_t)p;
    }
}

bool syndrome_ok(const LdpcCode& code, const uint8_t* c) {
    for (int chk = 0; chk < code.m; ++chk) {
        int parity = 0;
        for (int e = code.chk_ptr[chk]; e < code.chk_ptr[chk + 1]; ++e) parity ^= c[code.edge_var[e]];
        if (parity) return false;
    }
    return true;
}