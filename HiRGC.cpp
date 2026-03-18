#include <iostream>
#include <fstream>
#include <vector>
#include <string>

struct FastaData {
    std::string id;
    std::vector<std::string> sequences;
};

FastaData read_fasta(const std::string& filename) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    FastaData data;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            data.id = line.substr(1);
        } else {
            data.sequences.push_back(line);
        }
    }

    file.close();
    return data;
}

int main() {
    try {
        FastaData fasta = read_fasta("genomic.fna");

        std::cout << "ID: " << fasta.id << std::endl;
        std::cout << "Broj linija sekvence: " << fasta.sequences.size() << std::endl;

    } catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
    }

    return 0;
}