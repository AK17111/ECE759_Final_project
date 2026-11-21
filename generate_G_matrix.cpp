#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm> 

using namespace std;

// Function to parse an ALIST file into a dense 2D H-matrix
std::vector<std::vector<int>> parseAlistToHMatrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open ALIST file: " << filename << "\n";
        return {};
    }

    std::string line;
    int num_cols = 0; 
    int num_rows = 0;

    // 1. Read Line 1: Dimensions (V-nodes/Columns and C-nodes/Rows)
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        ss >> num_cols >> num_rows;
    }

    // 2. Initialize the 2D H matrix with all zeros
    std::vector<std::vector<int>> H(num_rows, std::vector<int>(num_cols, 0));

    // 3. Skip Lines 2 through 4 (Max weights, column weights, row weights)
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }

    // 4. Read the Variable Node block (1 line per column)
    for (int col = 0; col < num_cols; ++col) {
        if (!std::getline(file, line)) {
            std::cerr << "Error: Reached end of file prematurely at column " << col << "\n";
            break;
        }

        std::stringstream ss(line);
        int connected_row;
        
        while (ss >> connected_row) {
            if (connected_row > 0) { 
                H[connected_row - 1][col] = 1;
            }
        }
    }

    file.close();
    return H;
}

// Function to perform GF(2) Gauss-Jordan Elimination to derive G
vector<vector<int>> generateGMatrix(vector<vector<int>> H) {
    int M = H.size();       
    int N = H[0].size();    
    int K = N - M;          

    cout << "Starting GF(2) Gauss-Jordan Elimination...\n";

    // Phase 1: Transform H into Standard Form [P^T | I_M]
    for (int i = 0; i < M; ++i) {
        int target_col = K + i;
        int pivot_row = -1;

        for (int j = i; j < M; ++j) {
            if (H[j][target_col] == 1) {
                pivot_row = j;
                break;
            }
        }

        if (pivot_row == -1) {
            cerr << "Error: Matrix requires column permutations at column " << target_col << ". Not in standard form.\n";
            return {};
        }

        if (pivot_row != i) {
            swap(H[i], H[pivot_row]);
        }

        for (int j = 0; j < M; ++j) {
            if (j != i && H[j][target_col] == 1) {
                for (int col = 0; col < N; ++col) {
                    H[j][col] ^= H[i][col]; 
                }
            }
        }
    }

    // Phase 2: Extract P^T and Transpose it to P
    vector<vector<int>> P(K, vector<int>(M, 0));
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < K; ++j) {
            P[j][i] = H[i][j]; 
        }
    }

    // Phase 3: Construct the Systematic Generator Matrix G = [I_K | P]
    vector<vector<int>> G(K, vector<int>(N, 0));
    for (int i = 0; i < K; ++i) {
        G[i][i] = 1; 
        for (int j = 0; j < M; ++j) {
            G[i][K + j] = P[i][j];
        }
    }

    return G;
}

// Function to save a 2D matrix to a text file
void saveMatrixToFile(const string& filename, const vector<vector<int>>& mat) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error: Could not create output file: " << filename << "\n";
        return;
    }

    for (int i = 0; i < mat.size(); ++i) {
        for (int j = 0; j < mat[0].size(); ++j) {
            outfile << mat[i][j] << (j == mat[0].size() - 1 ? "" : " ");
        }
        outfile << "\n";
    }

    outfile.close();
    cout << "Successfully saved " << mat.size() << "x" << mat[0].size() << " matrix to " << filename << "\n";
}

int main() {
    string input_filename = "96.33.964.txt"; 
    string output_filename = "G_matrix_output.txt";
    
    cout << "Reading H Matrix from " << input_filename << "...\n";
    vector<vector<int>> H = parseAlistToHMatrix(input_filename);

    if (H.empty()) {
        return 1; 
    }
    cout << "Successfully generated " << H.size() << "x" << H[0].size() << " H-Matrix.\n\n";

    vector<vector<int>> G = generateGMatrix(H);

    if (!G.empty()) {
        saveMatrixToFile(output_filename, G);
    }

    return 0;
}
