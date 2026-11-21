#include <iostream>
#include <fstream>
#include <vector>
#include <random>

using namespace std;

// Function to read the 1D input array
vector<int> readInputData(const string& filename, int size) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << " ";
        return {};
    }
    vector<int> data(size);
    for (int i = 0; i < size; ++i) {
        file >> data[i];
    }
    file.close();
    return data;
}

// Function to read the 2D Generator Matrix
vector<vector<int>> readGMatrix(const string& filename, int rows, int cols) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << " ";
        return {};
    }
    vector<vector<int>> G(rows, vector<int>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file >> G[i][j];
        }
    }
    file.close();
    return G;
}

int main() {
    int K = 63; // Raw data bits
    int N = 96; // Total codeword bits
    
    // 1. Load Data
    cout << "Loading data...";
    vector<int> input_data = readInputData("random_input_63bit.txt", K);
    vector<vector<int>> G = readGMatrix("G_matrix_output.txt", K, N);
    
    if (input_data.empty() || G.empty()) {
        cerr << "Failed to load input files. Ensure they exist.";
        return 1;
    }

    // 2. Encode: Matrix Multiplication in GF(2) (Codeword = Input * G)
    cout << "Encoding payload...";
    vector<int> codeword(N, 0);
    for (int j = 0; j < N; ++j) {
        int bit_sum = 0;
        for (int i = 0; i < K; ++i) {
            // GF(2) multiplication is bitwise AND, addition is bitwise XOR
            bit_sum ^= (input_data[i] & G[i][j]); 
        }
        codeword[j] = bit_sum;
    }

    // 3. Modulate (BPSK) and Add AWGN
    cout << "Simulating AWGN Channel...";
    double noise_std_dev = 0.8; // Adjust this to change the Signal-to-Noise Ratio
    double variance = noise_std_dev * noise_std_dev;
    
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> awgn(0.0, noise_std_dev);

    vector<double> noisy_llrs(N);
    for (int i = 0; i < N; ++i) {
        // BPSK Mapping: 0 -> +1.0, 1 -> -1.0
        double voltage = (codeword[i] == 0) ? 1.0 : -1.0;
        
        // Add random Gaussian noise
        double received_signal = voltage + awgn(gen);
        
        // Calculate Log-Likelihood Ratio (LLR)
        noisy_llrs[i] = (2.0 * received_signal) / variance;
    }

    // 4. Save LLRs to text file for the decoder
    string output_filename = "noisy_llrs_input.txt";
    ofstream outfile(output_filename);
    if (outfile.is_open()) {
        for (int i = 0; i < N; ++i) {
            outfile << noisy_llrs[i] << (i == N - 1 ? "" : " ");
        }
        outfile << " ";
        outfile.close();
        cout << "Successfully saved " << N << " noisy LLRs to " << output_filename << " ";
    } else {
        cerr << "Error: Could not write to " << output_filename << " ";
    }

    return 0;
}
