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

struct PreprocessedData {
    string L1;
    string L2;
    string L3;

    vector<int> low_pos;
    vector<int> low_len;

    vector<int> N_pos;
    vector<int> N_len;

    vector<int> oth_pos;
    vector<char> oth_ch;
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

PreprocessedData preprocess_sequence(const string& input) {
    PreprocessedData result;
    bool in_lower = false;
    int start_lower = 0;

    for (int i = 0; i < input.size(); i++) {
        char c = input[i];

        if (islower(static_cast<unsigned char>(c))) {
            if (!in_lower) {
                in_lower = true;
                start_lower = i;
            }
        } else {
            if (in_lower) {
                result.low_pos.push_back(start_lower);
                result.low_len.push_back(i - start_lower);
                in_lower = false;
            }
        }

        result.L1 += toupper(static_cast<unsigned char>(c));
    }

    if (in_lower) {
        result.low_pos.push_back(start_lower);
        result.low_len.push_back(input.size() - start_lower);
    }

    bool in_N = false;
    int start_N = 0;

    for (int i = 0; i < result.L1.size(); i++) {
        char c = result.L1[i];

        if (c == 'N') {
            if (!in_N) {
                in_N = true;
                start_N = i;
            }
        } else {
            if (in_N) {
                result.N_pos.push_back(start_N);
                result.N_len.push_back(i - start_N);
                in_N = false;
            }
            result.L2 += c;
        }
    }

    if (in_N) {
        result.N_pos.push_back(start_N);
        result.N_len.push_back(result.L1.size() - start_N);
    }

    for (int i = 0; i < result.L2.size(); i++) {
        char upper = result.L2[i];

        if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') {
            result.L3 += upper;
        }
        else {
            result.oth_pos.push_back(i);
            result.oth_ch.push_back(upper);
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

        PreprocessedData clean_genome = preprocess_sequence(genome);
        cout << "L1 duljina: " << clean_genome.L1.size() << endl;
        cout << "L2 duljina: " << clean_genome.L2.size() << endl;
        cout << "L3 duljina: " << clean_genome.L3.size() << endl;

        cout << "Broj lower-case intervala: " << clean_genome.low_pos.size() << endl;
        cout << "Broj N intervala: " << clean_genome.N_pos.size() << endl;
        cout << "Broj ostalih znakova: " << clean_genome.oth_pos.size() << endl;

        string encoded = to_binary(clean_genome.L3);
        //cout << "Encoded (2-bit string): " << encoded << endl;
        cout << "Broj bitova: " << encoded.size() << endl;

    } catch (const exception& e) {
        cout << e.what() << endl;
    }

    return 0;
}