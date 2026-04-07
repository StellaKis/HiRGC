#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdint>
#include <cstdlib>

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

struct KTupleData {
    vector<int> values;
    vector<int> positions;
};

struct HashTable {
    vector<int> h;      
    vector<int> p;      
};

struct Match {
    int pos;    
    int len;    
};

struct Mismatch {
    int start;
    string seq;
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

vector<int> to_binary(const string& sequence) {
    vector<int> encoded;

    for (char base : sequence) {
        if (base == 'A') encoded.push_back(0);
        else if (base == 'C') encoded.push_back(1);
        else if (base == 'G') encoded.push_back(2);
        else if (base == 'T') encoded.push_back(3);
        else throw runtime_error("Invalid base");
    }

    return encoded;
}

char decode_base(int x) {
    if (x == 0) return 'A';
    if (x == 1) return 'C';
    if (x == 2) return 'G';
    if (x == 3) return 'T';
    throw runtime_error("Invalid base");
}

string decode_sequence(const vector<int>& encoded, int start, int end) {
    string result;

    for (int i = start; i <= end; i++) {
        result += decode_base(encoded[i]);
    }

    return result;
}

KTupleData generate_k_tuples(const vector<int>& encoded, int k) {
    KTupleData data;

    if (k <= 0 || k > encoded.size()) {
        return data;
    }

    for (int i = 0; i <= encoded.size() - k; i++) {
        int value = 0;
        int power = 1;

        for (int j = 0; j < k; j++) {
            value += encoded[i + j] * power;
            power *= 4;
        }

        data.values.push_back(value);
        data.positions.push_back(i);
    }

    return data;
}

HashTable generate_hash_table(const vector<int>& values, int s) {
    HashTable table;

    int n = values.size();

    table.h.assign(s, -1);   // inicijalno prazno
    table.p.assign(n, -1);   // next pointer

    for (int i = 0; i < n; i++) {
        int hash = values[i] % s;

        table.p[i] = table.h[hash]; // p(i) = h()
        table.h[hash] = i;          // h() = i
    }

    return table;
}

void greedy_matching(
    const vector<int>& ref_encoded,
    const vector<int>& target_encoded,
    int k,
    int s,
    vector<Match>& matches,
    vector<Mismatch>& mismatches
) {
    matches.clear();
    mismatches.clear();

    int nt = target_encoded.size();

    KTupleData ref_tuples = generate_k_tuples(ref_encoded, k);
    HashTable table = generate_hash_table(ref_tuples.values, s);

    int i = 0;
    int p_star = 0;

    while (i <= nt - k) {
        int Vt = 0;
        int power = 1;

        for (int j = 0; j < k; j++) {
            Vt += target_encoded[i + j] * power;
            power *= 4;
        }

        int bucket = Vt % s;
        int j = table.h[bucket];

        int p_max = -1;
        int l_max = 0;

        while (j != -1) {
            int l = 0;
            int ref_pos = ref_tuples.positions[j];

            while ((i + l < nt) &&
                   (ref_pos + l < ref_encoded.size()) &&
                   target_encoded[i + l] == ref_encoded[ref_pos + l]) {
                l++;
            }

            if (l >= k && l > l_max) {
                p_max = ref_pos;
                l_max = l;
            }

            j = table.p[j];
        }

        if (l_max > 0) {
            if (i > p_star) {
                Mismatch mm;
                mm.start = p_star;
                mm.seq = decode_sequence(target_encoded, p_star, i - 1);
                mismatches.push_back(mm);
            }

            Match m;
            m.pos = p_max;
            m.len = l_max;
            matches.push_back(m);

            p_star = i + l_max;
        }

        i = i + l_max + 1;
    }

    if (p_star < nt) {
        Mismatch mm;
        mm.start = p_star;
        mm.seq = decode_sequence(target_encoded, p_star, nt - 1);
        mismatches.push_back(mm);
    }
}

vector<int> delta_encode(const vector<int>& values) {
    vector<int> result;

    if (values.empty()) return result;

    result.push_back(values[0]);

    for (int i = 1; i < (int)values.size(); i++) {
        result.push_back(values[i] - values[i - 1]);
    }

    return result;
}

void rle_encode(const vector<size_t>& input, vector<int>& values, vector<int>& counts) {
    values.clear();
    counts.clear();

    if (input.empty()) return;

    int current = (int)input[0];
    int count = 1;

    for (int i = 1; i < (int)input.size(); i++) {
        if ((int)input[i] == current) {
            count++;
        } else {
            values.push_back(current);
            counts.push_back(count);
            current = (int)input[i];
            count = 1;
        }
    }

    values.push_back(current);
    counts.push_back(count);
}

vector<char> build_other_alphabet(const vector<char>& oth_ch) {
    vector<char> alphabet;

    for (char c : oth_ch) {
        bool found = false;

        for (char a : alphabet) {
            if (a == c) {
                found = true;
                break;
            }
        }

        if (!found) {
            alphabet.push_back(c);
        }
    }

    return alphabet;
}

vector<int> encode_other_chars(const vector<char>& oth_ch, const vector<char>& alphabet) {
    vector<int> codes;

    for (char c : oth_ch) {
        bool found = false;

        for (int i = 0; i < (int)alphabet.size(); i++) {
            if (alphabet[i] == c) {
                codes.push_back(i);
                found = true;
                break;
            }
        }

        if (!found) {
            throw runtime_error("Other character not found in alphabet");
        }
    }

    return codes;
}

void write_all_outputs(
    const string& matches_file,
    const string& mismatches_file,
    const string& aux_file,
    const vector<Match>& matches,
    const vector<Mismatch>& mismatches,
    const FastaData& fasta,
    const PreprocessedData& prep
) {
    vector<int> match_pos;
    vector<int> match_len;

    for (const auto& m : matches) {
        match_pos.push_back(m.pos);
        match_len.push_back(m.len);
    }

    vector<int> match_pos_delta = delta_encode(match_pos);

    vector<int> low_pos_delta = delta_encode(prep.low_pos);
    vector<int> N_pos_delta = delta_encode(prep.N_pos);
    vector<int> oth_pos_delta = delta_encode(prep.oth_pos);

    vector<int> seq_len_values;
    vector<int> seq_len_counts;
    rle_encode(fasta.seq_len, seq_len_values, seq_len_counts);

    vector<char> oth_alphabet = build_other_alphabet(prep.oth_ch);
    vector<int> oth_codes = encode_other_chars(prep.oth_ch, oth_alphabet);

    {
        ofstream out(matches_file);
        if (!out.is_open()) throw runtime_error("Cannot open matches file");

        out << matches.size() << "\n";
        for (int i = 0; i < (int)matches.size(); i++) {
            out << match_pos_delta[i] << " " << match_len[i] << "\n";
        }
    }

    {
        ofstream out(mismatches_file);
        if (!out.is_open()) throw runtime_error("Cannot open mismatches file");

        out << mismatches.size() << "\n";
        for (const auto& mm : mismatches) {
            out << mm.start << " " << mm.seq << "\n";
        }
    }

    {
        ofstream out(aux_file);
        if (!out.is_open()) throw runtime_error("Cannot open auxiliary file");

        out << "ID\n" << fasta.id << "\n";

        out << "SEQ_LEN_RLE\n" << seq_len_values.size() << "\n";
        for (int i = 0; i < (int)seq_len_values.size(); i++) {
            out << seq_len_values[i] << " " << seq_len_counts[i] << "\n";
        }

        out << "LOWERCASE\n" << prep.low_pos.size() << "\n";
        for (int i = 0; i < (int)prep.low_pos.size(); i++) {
            out << low_pos_delta[i] << " " << prep.low_len[i] << "\n";
        }

        out << "N_INTERVALS\n" << prep.N_pos.size() << "\n";
        for (int i = 0; i < (int)prep.N_pos.size(); i++) {
            out << N_pos_delta[i] << " " << prep.N_len[i] << "\n";
        }

        out << "OTHER_ALPHABET\n" << oth_alphabet.size() << "\n";
        for (char c : oth_alphabet) {
            out << c << " ";
        }
        out << "\n";

        out << "OTHER_CHARS\n";
        out << prep.oth_pos.size() << "\n";
        for (int i = 0; i < (int)prep.oth_pos.size(); i++) {
            out << oth_pos_delta[i] << " " << oth_codes[i] << "\n";
        }
    }
}

int compress_with_7zip(const string& archive_name, const string& matches_file, const string& mismatches_file, const string& aux_file) {
    string command = "7z a -t7z -m0=PPMd -mx=9 \"" + archive_name + "\" \"" +
                     matches_file + "\" \"" +
                     mismatches_file + "\" \"" +
                     aux_file + "\"";

    return system(command.c_str());
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

        vector<int> target_encoded = to_binary(clean_genome.L3);
        //cout << "Encoded (2-bit string): " << encoded << endl;
        cout << "Broj bitova: " << target_encoded.size() << endl;

        // int k = 3;
        // KTupleData ktuples = generate_k_tuples(encoded, k);
        // cout << "Broj k-tupleova: " << ktuples.values.size() << endl;

        //for (int i = 0; i < ktuples.values.size(); i++) {
        //    cout << "Tuple value: " << ktuples.values[i]
        //        << " na poziciji " << ktuples.positions[i] << endl;
        //}

        // int s = 10; // probno
        // HashTable table = generate_hash_table(ktuples.values, s);

        // cout << "\nHeader (h):\n";
        // for (int i = 0; i < table.h.size(); i++) {
        //     cout << "h[" << i << "] = " << table.h[i] << endl;
        // }

        // cout << "\nBuckets (linked lists):\n";
        // for (int i = 0; i < table.h.size(); i++) {
        //     cout << "Bucket " << i << ": ";

        //     int idx = table.h[i];

        //     while (idx != -1) {
        //         cout << idx << " -> ";
        //         idx = table.p[idx];  // idi na sljedeći element
        //     }

        //     cout << "NULL" << endl;
        // }

        FastaData fasta_ref = read_fasta("genomic_ref.fna");

        string genome_ref = join_sequences(fasta_ref.sequences);

        PreprocessedData clean_genome_ref = preprocess_sequence(genome_ref);

        vector<int> ref_encoded = to_binary(clean_genome_ref.L3);
        //cout << "Encoded (2-bit string): " << encoded << endl;
        cout << "Broj bitova: " << ref_encoded.size() << endl;

        int k = 3;
        int s = 10;

        vector<Match> matches;
        vector<Mismatch> mismatches;

        greedy_matching(ref_encoded, target_encoded, k, s, matches, mismatches);

        cout << "Broj match zapisa: " << matches.size() << endl;
        cout << "Broj mismatch zapisa: " << mismatches.size() << endl;

        // for (const auto& m : matches) {
        //     cout << "Match: pos=" << m.pos << ", len=" << m.len << endl;
        // }

        // for (const auto& mm : mismatches) {
        //     cout << "Mismatch: [" << mm.start << ", " << mm.end << "]" << endl;
        // }

        write_all_outputs("matches.txt", "mismatches.txt", "auxiliary.txt", matches, mismatches, fasta, clean_genome);

        int zip_status = compress_with_7zip("compressed_output.7z", "matches.txt", "mismatches.txt", "auxiliary.txt");

        if (zip_status != 0) {
            throw runtime_error("7-Zip compression failed");
        }

        cout << "7-Zip PPMd compression completed." << endl;

    } catch (const exception& e) {
        cout << e.what() << endl;
    }

    return 0;
}