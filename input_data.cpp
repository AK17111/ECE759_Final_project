#include <iostream>
#include <fstream>
#include <vector>
#include <random>

int main() {
    int raw_data_size = 63; // 96 total bits - 33 parity checks
    std::vector<int> random_input(raw_data_size);

    // 1. Set up a high-quality hardware random number generator
    std::random_device rd;  
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<> distrib(0, 1); // Only generate 0s and 1s

    // 2. Generate exactly 63 random bits
    for (int i = 0; i < raw_data_size; ++i) {
        random_input[i] = distrib(gen);
    }

    // 3. Save the random bits to a text file
    std::string filename = "random_input_63bit.txt";
    std::ofstream outfile(filename);
    
    if (outfile.is_open()) {
        for (int i = 0; i < raw_data_size; ++i) {
            outfile << random_input[i] << " ";
        }
        outfile.close();
        std::cout << "Successfully saved 63 random bits to " << filename << "\n";
    } else {
        std::cerr << "Error: Could not create the file.\n";
    }

    return 0;
}
