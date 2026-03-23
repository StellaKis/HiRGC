#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdint>

using namespace std;

struct FastaData {
    string id;
    vector<string> sequences;
    vector<size_t> seq_len;
};

FastaData read_fasta(const string& filename) {
    ifstream file(filename);
    
    if (!file.is_open()) {
        throw runtime_error("Unable to open file");
    }

    FastaData data;
    string line;

    while (getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            data.id = line.substr(1);
        } else {
            data.sequences.push_back(line);
            data.seq_len.push_back(line.size());
        }
    }

    file.close();
    return data;
}

string join_sequences(const vector<string>& sequences) {
    string result;

    for (const auto& seq : sequences) {
        result += seq;
    }

    return result;
}

string preprocess_sequence(const string& input) {
    string result;
    result.reserve(input.size());

    for (char c : input) {
        char upper = toupper(c);

        if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') {
            result += upper;
        }
        else if (upper == 'N') {
            continue;
        }
    }

    return result;
}

string to_binary(const string& sequence) {
    string encoded = "";

    for (char base : sequence) {
        if (base == 'A') encoded += "00";
        else if (base == 'C') encoded += "01";
        else if (base == 'G') encoded += "10";
        else if (base == 'T') encoded += "11";
        else throw runtime_error("Invalid base");
    }

    return encoded;
}

int main() {
    try {
        FastaData fasta = read_fasta("genomic.fna");

        cout << "ID: " << fasta.id << endl;
        cout << "Broj linija sekvence: " << fasta.sequences.size() << endl;

        // if (!fasta.seq_len.empty()) {
        //     cout << "Duljina prve linije: " << fasta.seq_len[0] << endl;
        // } else {
        //     cout << "Sekvenca nema linija!" << endl;
        // }

        string genome = join_sequences(fasta.sequences);
        cout << "Duljina sekvence: " << genome.size() << endl;

        string clean_genome = preprocess_sequence(genome);
        cout << "Prociscena duljina: " << clean_genome.size() << endl;

        string encoded = to_binary(clean_genome);
        //cout << "Encoded (2-bit string): " << encoded << endl;
        cout << "Broj bitova: " << encoded.size() << endl;

    } catch (const exception& e) {
        cout << e.what() << endl;
    }

    return 0;
}