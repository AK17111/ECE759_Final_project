#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm> // for std::min

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
    // Dimensions: num_rows x num_cols
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
        
        // Extract each space-separated integer on this line
        while (ss >> connected_row) {
            // Some ALIST formats pad short lines with 0s; valid indices are always >= 1
            if (connected_row > 0) { 
                // CRITICAL: Subtract 1 to shift ALIST 1-based index to C++ 0-based index
                H[connected_row - 1][col] = 1;
            }
        }
    }

    // File stream closes automatically when it goes out of scope, but doing it explicitly is good practice
    file.close();
    
    return H;
}

int main() {
    // Specify your file name here
    std::string filename = "96.33.964.txt"; 
    
    // Generate the matrix
    std::vector<std::vector<int>> H = parseAlistToHMatrix(filename);

    if (H.empty()) {
        return 1; // Exit if parsing failed
    }

    // Print summary
    std::cout << "Successfully generated " << H.size() << "x" << H[0].size() << " H-Matrix.\n\n";

    // Print the top-left 15x15 block of the matrix to visually verify the 1s and 0s
    std::cout << "Top-left 15x15 block visualization:\n";
    for (int i = 0; i < (int)H.size(); ++i) {
        for (int j = 0; j < (int)H[0].size(); ++j) {
            std::cout << H[i][j] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
