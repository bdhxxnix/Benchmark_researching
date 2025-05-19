#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

void readFilesInDirectory(const std::string& directoryPath) {
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.path().extension() == ".tsv") { // Process only .tsv files
            std::ifstream file(entry.path());
            if (!file) {
                std::cerr << "Error opening file: " << entry.path() << std::endl;
                continue;
            }

            std::cout << "Reading file: " << entry.path() << std::endl;
            std::string line;
            while (std::getline(file, line)) {
                std::cout << line << std::endl; // Process file content
            }
            file.close();
        }
    }
}

int main() {
    std::string folder = "./data"; // Directory containing TSV files
    readFilesInDirectory(folder);
    return 0;
}
