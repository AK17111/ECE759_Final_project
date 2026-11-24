#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

// =============================================================================
// ECE 759 Final Project: LDPC Min-Sum Iterative Decoder
// =============================================================================
// This program completes the receiver side of the LDPC flow:
// 1. Parses the ALIST file (96.33.964.txt) to construct the Tanner Graph.
// 2. Reads the noisy channel LLRs from 'noisy_llrs_input.txt'.
// 3. Executes the Normalized Min-Sum (NMS) message-passing decoding algorithm.
// 4. Performs syndrome validation (s = H * c^T mod 2) for early stopping.
// 5. Extracts the original systematic payload and verifies error correction.
// =============================================================================

// Structure to store Tanner Graph connectivity
struct TannerGraph {
    int N = 0; // Number of Variable Nodes (codeword length)
    int M = 0; // Number of Check Nodes (parity equations)
    int K = 0; // Number of Information/Payload bits (N - M)

    // Check-to-Variable connections: chk_neighbors[i] = list of variable nodes connected to check i
    vector<vector<int>> chk_neighbors;

    // Variable-to-Check connections: var_neighbors[j] = list of check nodes connected to variable j
    vector<vector<int>> var_neighbors;

    // Dense parity-check matrix H (for syndrome verification)
    vector<vector<int>> H;
};

// -----------------------------------------------------------------------------
// Step 1: Parse the ALIST file to extract Tanner Graph structure
// -----------------------------------------------------------------------------
bool parseAlist(const string& filename, TannerGraph& graph) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open ALIST file: " << filename << "\n";
        return false;
    }

    string line;

    // Line 1 contains: <num_variable_nodes (N)> <num_check_nodes (M)>
    if (getline(file, line)) {
        stringstream ss(line);
        ss >> graph.N >> graph.M;
    } else {
        cerr << "Error: Empty ALIST file.\n";
        return false;
    }

    graph.K = graph.N - graph.M;

    // Skip Lines 2, 3, 4 (max degrees and degree lists)
    for (int i = 0; i < 3; ++i) {
        getline(file, line);
    }

    // Initialize data structures
    graph.chk_neighbors.assign(graph.M, vector<int>());
    graph.var_neighbors.assign(graph.N, vector<int>());
    graph.H.assign(graph.M, vector<int>(graph.N, 0));

    // Next N lines: Check nodes connected to each Variable Node (1 line per column)
    for (int col = 0; col < graph.N; ++col) {
        if (!getline(file, line)) {
            cerr << "Error: Unexpected end of file while reading column " << col << "\n";
            return false;
        }

        stringstream ss(line);
        int connected_row;
        while (ss >> connected_row) {
            if (connected_row > 0) { // ALIST uses 1-based indexing; 0 means padding
                int row = connected_row - 1; // Convert to 0-based index
                graph.var_neighbors[col].push_back(row);
                graph.chk_neighbors[row].push_back(col);
                graph.H[row][col] = 1;
            }
        }
    }

    file.close();
    return true;
}

// -----------------------------------------------------------------------------
// Step 2: Read channel LLRs from 'noisy_llrs_input.txt'
// -----------------------------------------------------------------------------
vector<double> readChannelLLRs(const string& filename, int size) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << "\n";
        return {};
    }

    vector<double> llrs(size);
    for (int i = 0; i < size; ++i) {
        if (!(file >> llrs[i])) {
            cerr << "Error: Not enough LLR values in " << filename << "\n";
            return {};
        }
    }
    file.close();
    return llrs;
}

// -----------------------------------------------------------------------------
// Step 3: Read original raw input payload to verify bit correction
// -----------------------------------------------------------------------------
vector<int> readOriginalPayload(const string& filename, int size) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Warning: Could not open " << filename << " (Payload verification skipped).\n";
        return {};
    }

    vector<int> payload(size);
    for (int i = 0; i < size; ++i) {
        if (!(file >> payload[i])) {
            cerr << "Warning: Could not read " << size << " bits from " << filename << "\n";
            return {};
        }
    }
    file.close();
    return payload;
}

// -----------------------------------------------------------------------------
// Helper: Sign function for floating point LLRs (+1 for >= 0, -1 for < 0)
// -----------------------------------------------------------------------------
inline int sign(double val) {
    return (val >= 0.0) ? 1 : -1;
}

// -----------------------------------------------------------------------------
// Step 4: Iterative Min-Sum Decoding Algorithm
// -----------------------------------------------------------------------------
struct DecodeResult {
    bool converged;
    int iterations_taken;
    vector<int> decoded_codeword;
    vector<int> decoded_payload;
};

DecodeResult runMinSumDecoder(const TannerGraph& graph, 
                             const vector<double>& channel_llrs, 
                             int max_iterations = 50, 
                             double alpha = 0.8) {
    int N = graph.N;
    int M = graph.M;
    int K = graph.K;

    // Messages on Tanner graph edges:
    // L_v2c[i][j]: Message from Variable Node j to Check Node i
    // L_c2v[i][j]: Message from Check Node i to Variable Node j
    vector<vector<double>> L_v2c(M, vector<double>(N, 0.0));
    vector<vector<double>> L_c2v(M, vector<double>(N, 0.0));

    // -------------------------------------------------------------------------
    // Phase 4.1: Initialization (Iteration 0)
    // Every variable node sends its initial channel LLR to all connected checks
    // -------------------------------------------------------------------------
    for (int j = 0; j < N; ++j) {
        for (int i : graph.var_neighbors[j]) {
            L_v2c[i][j] = channel_llrs[j];
        }
    }

    vector<double> L_total(N, 0.0);
    vector<int> c_hat(N, 0);
    bool converged = false;
    int iter = 0;

    cout << "\n============================================================\n";
    cout << "             Starting Min-Sum Decoding Loop                 \n";
    cout << "============================================================\n";
    cout << "Max Iterations : " << max_iterations << "\n";
    cout << "Min-Sum Scale Factor (alpha): " << alpha << "\n\n";

    for (iter = 1; iter <= max_iterations; ++iter) {
        // ---------------------------------------------------------------------
        // Phase 4.2: Check Node Update (Horizontal Step)
        // Each check node calculates parity consistency for all its neighbors.
        // The message to neighbor j must EXCLUDE j's own incoming message (Extrinsic).
        // ---------------------------------------------------------------------
        for (int i = 0; i < M; ++i) {
            const auto& neighbors = graph.chk_neighbors[i];
            int degree = neighbors.size();

            for (int k = 0; k < degree; ++k) {
                int target_var = neighbors[k];

                int prod_sign = 1;
                double min_magnitude = 1e9; // Large positive number

                // Scan all OTHER connected variable nodes
                for (int m = 0; m < degree; ++m) {
                    if (m == k) continue; // Exclude target variable node!

                    int other_var = neighbors[m];
                    double val = L_v2c[i][other_var];

                    prod_sign *= sign(val);
                    min_magnitude = min(min_magnitude, fabs(val));
                }

                // Check-to-Variable Message:
                // Normalized Min-Sum applies scaling factor alpha (e.g. 0.8) to counter overconfidence
                L_c2v[i][target_var] = alpha * prod_sign * min_magnitude;
            }
        }

        // ---------------------------------------------------------------------
        // Phase 4.3: Variable Node Update (Vertical Step) & Total Accumulation
        // Each variable node aggregates all check messages + original channel LLR.
        // ---------------------------------------------------------------------
        for (int j = 0; j < N; ++j) {
            // First, calculate the total sum of all incoming check messages
            double check_messages_sum = 0.0;
            for (int i : graph.var_neighbors[j]) {
                check_messages_sum += L_c2v[i][j];
            }

            // Total accumulated belief (combines channel observation + all check advice)
            L_total[j] = channel_llrs[j] + check_messages_sum;

            // Extrinsic message to check node i: exclude check i's own advice
            for (int i : graph.var_neighbors[j]) {
                L_v2c[i][j] = L_total[j] - L_c2v[i][j];
            }

            // -----------------------------------------------------------------
            // Phase 4.4: Hard Decision Slicing
            // LLR >= 0.0 -> Bit 0 (positive voltage favored)
            // LLR <  0.0 -> Bit 1 (negative voltage favored)
            // -----------------------------------------------------------------
            c_hat[j] = (L_total[j] >= 0.0) ? 0 : 1;
        }

        // ---------------------------------------------------------------------
        // Phase 4.5: Syndrome Parity Check (s = H * c_hat^T mod 2)
        // Check if every parity equation is satisfied.
        // ---------------------------------------------------------------------
        int unsatisfied_checks = 0;
        for (int i = 0; i < M; ++i) {
            int parity_sum = 0;
            for (int var_idx : graph.chk_neighbors[i]) {
                parity_sum ^= c_hat[var_idx];
            }
            if (parity_sum != 0) {
                unsatisfied_checks++;
            }
        }

        // Print progress
        cout << "Iteration " << setw(2) << iter 
             << " | Unsatisfied Parity Checks: " << setw(2) << unsatisfied_checks;

        if (unsatisfied_checks == 0) {
            cout << "  ===> [CONVERGED! Halting Early]\n";
            converged = true;
            break;
        } else {
            cout << "\n";
        }
    }

    if (!converged) {
        cout << "\n[!] Decoder reached maximum iterations without satisfying all parity checks.\n";
    }

    // Extract systematic payload (the first K bits)
    vector<int> payload(c_hat.begin(), c_hat.begin() + K);

    return {converged, iter, c_hat, payload};
}

// -----------------------------------------------------------------------------
// Main Function
// -----------------------------------------------------------------------------
int main() {
    cout << "============================================================\n";
    cout << "             LDPC Pipeline: Receiver Decoder                \n";
    cout << "============================================================\n";

    // 1. Parse Tanner graph from the ALIST file
    string alist_file = "96.33.964.txt";
    cout << "Loading Tanner Graph from " << alist_file << "... ";
    TannerGraph graph;
    if (!parseAlist(alist_file, graph)) {
        return 1;
    }
    cout << "Done!\n";
    cout << "  - Codeword length (N): " << graph.N << " bits\n";
    cout << "  - Parity checks   (M): " << graph.M << " equations\n";
    cout << "  - Systematic data (K): " << graph.K << " bits\n";

    // 2. Read noisy channel LLRs generated by val.cpp
    string llr_file = "noisy_llrs_input.txt";
    cout << "Loading Channel LLRs from " << llr_file << "... ";
    vector<double> channel_llrs = readChannelLLRs(llr_file, graph.N);
    if (channel_llrs.empty()) {
        return 1;
    }
    cout << "Done! (" << channel_llrs.size() << " LLRs loaded)\n";

    // Analyze raw channel errors prior to decoding (uncorrected bit flips)
    vector<int> raw_channel_decisions(graph.N);
    for (int j = 0; j < graph.N; ++j) {
        raw_channel_decisions[j] = (channel_llrs[j] >= 0.0) ? 0 : 1;
    }

    // 3. Load original raw payload to measure error correction success
    string input_file = "random_input_63bit.txt";
    vector<int> orig_payload = readOriginalPayload(input_file, graph.K);

    int initial_payload_errors = 0;
    if (!orig_payload.empty()) {
        for (int i = 0; i < graph.K; ++i) {
            if (raw_channel_decisions[i] != orig_payload[i]) {
                initial_payload_errors++;
            }
        }
        cout << "Initial channel errors in payload (before decoding): " 
             << initial_payload_errors << " / " << graph.K << " bits\n";
    }

    // 4. Run the Min-Sum Iterative Decoder
    DecodeResult result = runMinSumDecoder(graph, channel_llrs, 50, 0.8);

    // 5. Verify Results
    cout << "\n============================================================\n";
    cout << "                     Decoding Summary                       \n";
    cout << "============================================================\n";
    cout << "Status: " << (result.converged ? "SUCCESS (Valid Codeword Found)" : "FAILURE") << "\n";
    cout << "Total Iterations Taken: " << result.iterations_taken << "\n";

    if (!orig_payload.empty()) {
        int final_errors = 0;
        for (int i = 0; i < graph.K; ++i) {
            if (result.decoded_payload[i] != orig_payload[i]) {
                final_errors++;
            }
        }
        cout << "Payload Bit Errors Remaining: " << final_errors << " / " << graph.K << "\n";
        if (final_errors == 0) {
            cout << "===> ALL ERRORS SUCCESSFULLY CORRECTED! Payload is 100% restored.\n";
        } else {
            cout << "===> Some bit errors remain.\n";
        }
    }

    // 6. Save Decoded Codeword and Payload to text files
    string codeword_out = "decoded_codeword.txt";
    ofstream f_cw(codeword_out);
    if (f_cw.is_open()) {
        for (int bit : result.decoded_codeword) f_cw << bit << " ";
        f_cw.close();
        cout << "\nSaved decoded codeword (" << graph.N << " bits) to " << codeword_out << "\n";
    }

    string payload_out = "decoded_payload.txt";
    ofstream f_pl(payload_out);
    if (f_pl.is_open()) {
        for (int bit : result.decoded_payload) f_pl << bit << " ";
        f_pl.close();
        cout << "Saved decoded payload (" << graph.K << " bits) to " << payload_out << "\n";
    }

    cout << "============================================================\n";
    return 0;
}
