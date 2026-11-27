//   ./ldpc test                                                        run all correctness tests
//   ./ldpc gencode   --n 2048 --out codes/peg_2048_3_6.alist           make an H matrix (PEG)
//   ./ldpc genframes --code codes/peg_2048_3_6.alist --ebno 1.75 --frames 1000 --out data/f.bin
//   ./ldpc sweep     --code codes/peg_2048_3_6.alist --sweep 1.0:0.25:2.5 --frames 2000 [--csv results/x.csv]
//   ./ldpc decode    --code codes/peg_2048_3_6.alist --in data/f.bin [--out results/serial.dec]
//   ./ldpc compare   results/serial.dec results/omp.dec
//
// Decoder options (sweep / decode): --iters 20  --alpha 0.75  --no-early-stop  --repeat 3  --seed 7
// OpenMP options   (sweep / decode): --threads 4  --schedule dynamic|static   (no --threads = serial)
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>

#include "ldpc.hpp"

// measuring
struct Stats {
    long info_bit_errors = 0, frame_errors = 0, total_iters = 0;
    double seconds = 0;
    int threads = 0;   // 0 = serial
};

// Decode all codewords (timed), then count errors against the original bits.
static Stats run_and_measure(const LdpcCode& code, const FrameSet& fs, const DecodeParams& p,
                             int repeat, std::vector<uint8_t>& bits, std::vector<int>& iters) {
    bits.assign((size_t)fs.frames * code.n, 0);
    iters.assign(fs.frames, 0);
    Stats s;
    s.seconds = 1e30;
    s.threads = p.threads;
    for (int r = 0; r < repeat; ++r) {  // keep the fastest run
        auto t0 = std::chrono::steady_clock::now();
        if (p.threads > 0)
            decode_batch_omp(code, fs.frames, fs.llr.data(), bits.data(), iters.data(), p);
        else
            decode_batch_serial(code, fs.frames, fs.llr.data(), bits.data(), iters.data(), p);
        auto t1 = std::chrono::steady_clock::now();
        s.seconds = std::min(s.seconds, std::chrono::duration<double>(t1 - t0).count());
    }
    for (int f = 0; f < fs.frames; ++f) {
        long errs = 0;
        for (int v : code.info_pos) errs += bits[(size_t)f * code.n + v] != fs.truth[(size_t)f * code.n + v];
        s.info_bit_errors += errs;
        s.frame_errors += errs > 0;
        s.total_iters += iters[f];
    }
    return s;
}

static void print_header() {
    std::printf("%8s %8s %12s %10s %9s %10s %12s %12s\n", "EbN0_dB", "frames", "BER", "FER",
                "avg_iter", "time_ms", "info_Mbps", "coded_Mbps");
}

// Print one result line and return the same numbers as a CSV line.
static std::string print_row(const LdpcCode& code, float ebno, int frames, const Stats& s) {
    double ber = (double)s.info_bit_errors / ((double)frames * code.k);
    double fer = (double)s.frame_errors / frames;
    double info_mbps = (double)frames * code.k / s.seconds / 1e6;
    double coded_mbps = (double)frames * code.n / s.seconds / 1e6;
    std::printf("%8.2f %8d %12.3e %10.3e %9.2f %10.2f %12.3f %12.3f\n", ebno, frames, ber, fer,
                (double)s.total_iters / frames, s.seconds * 1e3, info_mbps, coded_mbps);
    std::ostringstream o;
    o << code.n << "," << code.k << "," << ebno << "," << frames << "," << s.info_bit_errors << ","
      << s.frame_errors << "," << ber << "," << fer << "," << (double)s.total_iters / frames << ","
      << s.seconds * 1e3 << "," << info_mbps << "," << coded_mbps << "," << s.threads;
    return o.str();
}

static void append_csv(const std::string& path, const std::vector<std::string>& rows) {
    bool exists = (bool)std::ifstream(path);
    std::ofstream o(path, std::ios::app);
    if (!o) throw std::runtime_error("cannot write " + path + " (does the folder exist?)");
    if (!exists)
        o << "n,k,ebno_db,frames,info_bit_errors,frame_errors,ber,fer,avg_iters,time_ms,info_mbps,coded_mbps,threads\n";
    for (auto& r : rows) o << r << "\n";
}

static LdpcCode load_code(const Args& a) {
    LdpcCode code = load_alist(a.need("--code"));
    build_encoder(code);
    return code;
}

static DecodeParams read_params(const Args& a, const LdpcCode& code) {
    DecodeParams p;
    p.max_iters = (int)a.num("--iters", 20);
    p.alpha = (float)a.real("--alpha", 0.75);
    p.early_stop = !a.has("--no-early-stop");
    p.threads = (int)a.num("--threads", 0);                  // 0 = serial
    p.dynamic_sched = a.str("--schedule", "dynamic") != "static";
    std::printf("code: n=%d m=%d edges=%d k=%d rate=%.4f | iters=%d alpha=%.2f early_stop=%d | %s\n",
                code.n, code.m, code.num_edges(), code.k, code.rate(), p.max_iters, p.alpha,
                (int)p.early_stop,
                p.threads > 0 ? (std::string("OpenMP threads=") + std::to_string(p.threads) +
                                 (p.dynamic_sched ? " schedule=dynamic" : " schedule=static")).c_str()
                              : "serial");
    return p;
}

//  commands
static int cmd_gencode(const Args& a) {
    LdpcCode code = make_peg_code((int)a.num("--n", 96), (int)a.num("--dv", 3), (int)a.num("--dc", 6),
                                  (uint64_t)a.num("--seed", 1));
    build_encoder(code);
    save_alist(code, a.need("--out"));
    std::printf("n=%d m=%d edges=%d k=%d rate=%.4f 4-cycles=%ld -> %s\n", code.n, code.m,
                code.num_edges(), code.k, code.rate(), count_4cycles(code), a.need("--out").c_str());
    return 0;
}

static int cmd_genframes(const Args& a) {
    LdpcCode code = load_code(a);
    FrameSet fs = make_frames(code, (int)a.num("--frames", 1000), a.real("--ebno", 2.0),
                              (uint64_t)a.num("--seed", 7));
    write_frames(fs, a.need("--out"));
    std::printf("wrote %d codewords of %d bits at Eb/N0=%.2f dB -> %s\n", fs.frames, fs.n, fs.ebno_db,
                a.need("--out").c_str());
    return 0;
}

static int cmd_sweep(const Args& a) {
    LdpcCode code = load_code(a);
    DecodeParams p = read_params(a, code);
    double lo = 1.0, step = 0.5, hi = 3.0;
    std::sscanf(a.str("--sweep", "1.0:0.5:3.0").c_str(), "%lf:%lf:%lf", &lo, &step, &hi);
    int frames = (int)a.num("--frames", 1000), repeat = std::max(1, (int)a.num("--repeat", 1));
    uint64_t seed = (uint64_t)a.num("--seed", 7);
    std::vector<uint8_t> bits;
    std::vector<int> iters;
    std::vector<std::string> csv;
    print_header();
    int idx = 0;
    for (double eb = lo; eb <= hi + 1e-9; eb += step, ++idx) {
        FrameSet fs = make_frames(code, frames, eb, seed + idx);
        csv.push_back(print_row(code, (float)eb, frames, run_and_measure(code, fs, p, repeat, bits, iters)));
    }
    if (a.has("--csv")) append_csv(a.need("--csv"), csv);
    return 0;
}

static int cmd_decode(const Args& a) {
    LdpcCode code = load_code(a);
    DecodeParams p = read_params(a, code);
    FrameSet fs = read_frames(a.need("--in"));
    if (fs.n != code.n) throw std::runtime_error("frame file does not match code length");
    std::vector<uint8_t> bits;
    std::vector<int> iters;
    print_header();
    std::string row = print_row(code, fs.ebno_db, fs.frames,
                                run_and_measure(code, fs, p, std::max(1, (int)a.num("--repeat", 1)), bits, iters));
    if (a.has("--out")) write_decoded(a.need("--out"), code.n, fs.frames, bits, iters);
    if (a.has("--csv")) append_csv(a.need("--csv"), {row});
    return 0;
}

static int cmd_compare(int argc, char** argv) {
    if (argc != 4) throw std::runtime_error("usage: ./ldpc compare A.dec B.dec");
    int na, fa, nb, fb;
    std::vector<uint8_t> ba, bb;
    std::vector<int> ia, ib;
    read_decoded(argv[2], na, fa, ba, ia);
    read_decoded(argv[3], nb, fb, bb, ib);
    if (na != nb || fa != fb) { std::printf("shape mismatch\nDIFFERENT\n"); return 1; }
    long bit_diff = 0, frame_diff = 0, iter_diff = 0;
    for (int f = 0; f < fa; ++f) {
        long d = 0;
        for (int j = 0; j < na; ++j) d += ba[(size_t)f * na + j] != bb[(size_t)f * na + j];
        bit_diff += d;
        frame_diff += d > 0;
        iter_diff += ia[f] != ib[f];
    }
    std::printf("codewords=%d bits each=%d  differing bits=%ld  differing codewords=%ld  "
                "differing iteration counts=%ld\n", fa, na, bit_diff, frame_diff, iter_diff);
    bool same = bit_diff == 0 && iter_diff == 0;
    std::printf(same ? "IDENTICAL\n" : "DIFFERENT\n");
    return same ? 0 : 1;
}

// --- tests
static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("  %s  %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}

static int cmd_test() {
    std::printf("[1] 6-bit example: b2 misread as -1 gets fixed in one iteration\n");
    {
        LdpcCode code = code_from_rows(6, {{0, 1, 2}, {2, 3, 4}, {0, 4, 5}});
        float llr[6] = {4, 3, -1, 5, 2, 6};
        uint8_t bits[6];
        DecodeParams p;
        DecoderWorkspace ws(code);
        bool ok;
        int it = decode_frame(code, llr, bits, p, ws, &ok);
        bool zeros = true;
        for (uint8_t b : bits) zeros &= b == 0;
        check(ok && zeros, "decodes to all zeros");
        check(it == 1, "needs exactly 1 iteration");
        check(std::fabs(ws.c2v[2] - 2.25f) < 1e-6 && std::fabs(ws.c2v[3] - 1.5f) < 1e-6,
              "check messages to b2 are +2.25 and +1.5");
    }
    std::printf("[2] PDF Figure 3: inputs {+3.2,-1.1,-4.5} -> message +1.1 (alpha=1)\n");
    {
        LdpcCode code = code_from_rows(4, {{0, 1, 2, 3}});
        float llr[4] = {3.2f, -1.1f, -4.5f, 9.0f};
        uint8_t bits[4];
        DecodeParams p; p.alpha = 1.0f; p.max_iters = 1; p.early_stop = false;
        DecoderWorkspace ws(code);
        bool ok;
        decode_frame(code, llr, bits, p, ws, &ok);
        check(std::fabs(ws.c2v[3] - 1.1f) < 1e-6, "check message to v_j == +1.1");
    }
    std::printf("[3] PEG (3,6) code, n=96 (same shape as MacKay 96.33.964)\n");
    LdpcCode c96 = make_peg_code(96, 3, 6, 1);
    build_encoder(c96);
    {
        bool regular = true;
        for (int j = 0; j < c96.n; ++j) regular &= c96.var_ptr[j + 1] - c96.var_ptr[j] == 3;
        for (int c = 0; c < c96.m; ++c) regular &= c96.chk_ptr[c + 1] - c96.chk_ptr[c] == 6;
        check(c96.m == 48 && c96.num_edges() == 288 && regular, "48 checks, 288 edges, (3,6)-regular");
        check(count_4cycles(c96) == 0, "no 4-cycles (girth >= 6)");
        check(c96.k >= 48, "k >= 48 payload bits");
    }
    std::printf("[4] alist save/load round trip\n");
    {
        save_alist(c96, "/tmp/test96.alist");
        LdpcCode back = load_alist("/tmp/test96.alist");
        check(back.chk_ptr == c96.chk_ptr && back.edge_var == c96.edge_var, "identical graph");
    }
    std::printf("[5] encoder: 200 random payloads all satisfy H*c = 0\n");
    {
        LdpcCode c = make_peg_code(2048, 3, 6, 3);
        build_encoder(c);
        std::mt19937 rng(1);
        std::vector<uint8_t> info(c.k), cw(c.n);
        bool all = true;
        for (int t = 0; t < 200; ++t) {
            for (auto& b : info) b = rng() & 1;
            encode(c, info.data(), cw.data());
            all &= syndrome_ok(c, cw.data());
            for (int i = 0; i < c.k; ++i) all &= cw[c.info_pos[i]] == info[i];
        }
        check(all, "valid codewords, payload bits in place");

        std::printf("[6] noise-free codewords decode with zero changes and 0 iterations\n");
        DecodeParams p;
        DecoderWorkspace ws(c);
        std::vector<float> llr(c.n);
        std::vector<uint8_t> out(c.n);
        bool good = true;
        for (int t = 0; t < 20; ++t) {
            for (auto& b : info) b = rng() & 1;
            encode(c, info.data(), cw.data());
            for (int j = 0; j < c.n; ++j) llr[j] = cw[j] ? -10.0f : 10.0f;
            bool ok;
            int it = decode_frame(c, llr.data(), out.data(), p, ws, &ok);
            good &= ok && it == 0 && out == cw;
        }
        check(good, "20/20 clean codewords unchanged");

        std::printf("[7] codewords with 10 known wrong bits are fully corrected\n");
        int fixed = 0;
        for (int t = 0; t < 50; ++t) {
            for (auto& b : info) b = rng() & 1;
            encode(c, info.data(), cw.data());
            for (int j = 0; j < c.n; ++j) llr[j] = cw[j] ? -4.0f : 4.0f;
            for (int e = 0; e < 10; ++e) {
                int j = rng() % c.n;
                llr[j] = cw[j] ? 1.0f : -1.0f;  // wrong sign, low confidence
            }
            bool ok;
            decode_frame(c, llr.data(), out.data(), p, ws, &ok);
            fixed += ok && out == cw;
        }
        check(fixed == 50, "50/50 codewords corrected");
    }
    std::printf("[8] frame file round trip\n");
    {
        FrameSet fs = make_frames(c96, 5, 2.0, 11);
        write_frames(fs, "/tmp/test.bin");
        FrameSet back = read_frames("/tmp/test.bin");
        check(back.llr == fs.llr && back.truth == fs.truth, "identical frames");
    }
    std::printf(failures ? "\n%d TEST(S) FAILED\n" : "\nALL TESTS PASSED\n", failures);
    return failures ? 1 : 0;
}

// main
int main(int argc, char** argv) {
    try {
        std::string cmd = argc > 1 ? argv[1] : "";
        Args a(argc, argv);
        if (cmd == "test") return cmd_test();
        if (cmd == "gencode") return cmd_gencode(a);
        if (cmd == "genframes") return cmd_genframes(a);
        if (cmd == "sweep") return cmd_sweep(a);
        if (cmd == "decode") return cmd_decode(a);
        if (cmd == "compare") return cmd_compare(argc, argv);
        std::printf("usage: ./ldpc test | gencode | genframes | sweep | decode | compare  (see top of main.cpp)\n");
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}