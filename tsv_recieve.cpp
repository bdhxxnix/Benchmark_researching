#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>

int main() {
    std::ifstream infile("../data/UCRArchive_2018/CricketY/CricketY_TRAIN.tsv"); // Open the input TSV file
    if (!infile.is_open()) {
        std::cerr << "Error opening input file!\n";
        return 1;
    }

    std::ofstream outfile("CricketY_output.txt"); // Open the output file
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file!\n";
        return 1;
    }

    std::string line;

    // Read each line from the input file
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::vector<std::string> row;
        std::string value;

        // Read each column (split by tab '\t')
        while (std::getline(ss, value, '\t')) {
            row.push_back(value);
        }


        // Skip writing the first column (if there is at least one column)
        if (row.size() > 1) {
            for (size_t i = 1; i < row.size(); ++i) { // Start from index 1
                outfile << row[i];
                if (i < row.size() - 1) {
                    outfile << "\t"; // Add tab separator between columns
                }
                
            }
            outfile << "\n"; // New line for each row
        }
    }

    infile.close();
    outfile.close();

    
    return 0;
}
