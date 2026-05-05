#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <chrono>

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
    vector<uint64_t> values;
    vector<int> positions;
};

struct HashTable {
    vector<int> h;      
    vector<int> p;      
};

struct Match {
    int start_target; 
    int pos;    
    int len;    
};

struct Mismatch {
    int start;
    string seq;
};

struct AuxData {
    string id;

    vector<size_t> seq_len;

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
        uint64_t value = 0;
        uint64_t power = 1;

        for (int j = 0; j < k; j++) {
            value += encoded[i + j] * power;
            power *= 4;
        }

        data.values.push_back(value);
        data.positions.push_back(i);
    }

    return data;
}

HashTable generate_hash_table(const vector<uint64_t>& values, int s) {
    HashTable table;

    int n = (int)values.size();

    table.h.assign(s, -1);
    table.p.assign(n, -1);

    int mask = s - 1;

    for (int i = 0; i < n; i++) {
        int hash = (int)(values[i] & mask);

        table.p[i] = table.h[hash];
        table.h[hash] = i;
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

    int nr = (int)ref_encoded.size();
    int nt = (int)target_encoded.size();

    if (k <= 0 || nr < k || nt < k) {
        if (nt > 0) {
            Mismatch mm;
            mm.start = 0;
            mm.seq = decode_sequence(target_encoded, 0, nt - 1);
            mismatches.push_back(mm);
        }
        return;
    }

    KTupleData ref_tuples;
    int ref_tuple_count = nr - k + 1;

    uint64_t ref_value = 0;

    for (int j = k - 1; j >= 0; j--) {
        ref_value <<= 2;
        ref_value += ref_encoded[j];
    }

    ref_tuples.values.push_back(ref_value);
    ref_tuples.positions.push_back(0);

    int shift_bits = 2 * (k - 1);

    for (int i = 1; i < ref_tuple_count; i++) {
        ref_value >>= 2;
        ref_value += ((uint64_t)ref_encoded[i + k - 1] << shift_bits);

        ref_tuples.values.push_back(ref_value);
        ref_tuples.positions.push_back(i);
    }

    HashTable table = generate_hash_table(ref_tuples.values, s);

    int mask = s - 1;

    int i = 0;
    int p_star = 0;

    uint64_t Vt = 0;

    for (int j = k - 1; j >= 0; j--) {
        Vt <<= 2;
        Vt += target_encoded[j];
    }

    int current_hash_pos = 0;

    while (i <= nt - k) {
        while (current_hash_pos < i) {
            Vt >>= 2;
            Vt += ((uint64_t)target_encoded[current_hash_pos + k] << shift_bits);
            current_hash_pos++;
        }

        int bucket = (int)(Vt & mask);
        int j = table.h[bucket];

        int p_max = -1;
        int l_max = 0;

        while (j != -1) {
            if (ref_tuples.values[j] != Vt) {
                j = table.p[j];
                continue;
            }

            int ref_pos = ref_tuples.positions[j];

            int l = k;

            while ((i + l < nt) &&
                (ref_pos + l < nr) &&
                target_encoded[i + l] == ref_encoded[ref_pos + l]) {
                l++;
            }

            if (l > l_max) {
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
            m.start_target = i;
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
    vector<int> match_start;
    vector<int> match_pos;
    vector<int> match_len;

    for (const auto& m : matches) {
        match_start.push_back(m.start_target);
        match_pos.push_back(m.pos);
        match_len.push_back(m.len);
    }

    vector<int> match_start_delta = delta_encode(match_start);
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
            out << match_start_delta[i] << " " << match_pos_delta[i] << " " << match_len[i] << "\n";
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
    string command = "7z a -t7z -m0=PPMd -mx=5 \"" + archive_name + "\" \"" + matches_file + "\" \"" + mismatches_file + "\" \"" + aux_file + "\"";

    return system(command.c_str());
}

//dekompresija

int extract_with_7zip(const string& archive_name, const string& output_dir) {
    string command = "7z x \"" + archive_name + "\" -o\"" + output_dir + "\" -y";

    return system(command.c_str());
}

vector<int> delta_decode(const vector<int>& delta) {
    vector<int> result;

    if (delta.empty()) return result;

    result.push_back(delta[0]);

    for (int i = 1; i < (int)delta.size(); i++) {
        result.push_back(result[i - 1] + delta[i]);
    }

    return result;
}

vector<size_t> rle_decode(const vector<int>& values, const vector<int>& counts) {
    vector<size_t> result;

    for (int i = 0; i < (int)values.size(); i++) {
        for (int j = 0; j < counts[i]; j++) {
            result.push_back(values[i]);
        }
    }

    return result;
}

vector<char> decode_other_chars(const vector<int>& codes, const vector<char>& alphabet) {
    vector<char> result;

    for (int code : codes) {
        result.push_back(alphabet[code]);
    }

    return result;
}

vector<Match> read_matches(const string& filename) {
    ifstream in(filename);

    if (!in.is_open()) {
        throw runtime_error("Cannot open matches file");
    }

    int n;
    in >> n;

    vector<int> start_delta(n);
    vector<int> pos_delta(n);
    vector<int> len(n);

    for (int i = 0; i < n; i++) {
        in >> start_delta[i]
           >> pos_delta[i]
           >> len[i];
    }

    vector<int> start = delta_decode(start_delta);
    vector<int> pos = delta_decode(pos_delta);

    vector<Match> matches(n);

    for (int i = 0; i < n; i++) {
        matches[i].start_target = start[i];
        matches[i].pos = pos[i];
        matches[i].len = len[i];
    }

    return matches;
}

vector<Mismatch> read_mismatches(const string& filename) {
    ifstream in(filename);

    if (!in.is_open()) {
        throw runtime_error("Cannot open mismatches");
    }

    int n;
    in >> n;

    vector<Mismatch> mismatches(n);

    for (int i = 0; i < n; i++) {
        in >> mismatches[i].start
           >> mismatches[i].seq;
    }

    return mismatches;
}

AuxData read_auxiliary(const string& filename) {

    ifstream in(filename);

    if (!in.is_open()) {
        throw runtime_error("Cannot open auxiliary file");
    }

    AuxData aux;

    string token;

    // ID

    in >> token; // ID
    in.ignore();

    getline(in, aux.id);

    //SEQ_LEN_RLE

    in >> token;

    int rle_n;
    in >> rle_n;

    vector<int> values(rle_n);
    vector<int> counts(rle_n);

    for (int i = 0; i < rle_n; i++) {
        in >> values[i] >> counts[i];
    }

    aux.seq_len = rle_decode(values, counts);

    //LOWERCASE

    in >> token;

    int low_n;
    in >> low_n;

    vector<int> low_delta(low_n);

    aux.low_len.resize(low_n);

    for (int i = 0; i < low_n; i++) {
        in >> low_delta[i]
           >> aux.low_len[i];
    }

    aux.low_pos = delta_decode(low_delta);

    //N_INTERVALS

    in >> token;

    int N_n;
    in >> N_n;

    vector<int> N_delta(N_n);

    aux.N_len.resize(N_n);

    for (int i = 0; i < N_n; i++) {
        in >> N_delta[i]
           >> aux.N_len[i];
    }

    aux.N_pos = delta_decode(N_delta);

    //OTHER_ALPHABET

    in >> token;

    int alphabet_size;
    in >> alphabet_size;

    vector<char> alphabet(alphabet_size);

    for (int i = 0; i < alphabet_size; i++) {
        in >> alphabet[i];
    }

    //OTHER_CHARS

    in >> token;

    int oth_n;
    in >> oth_n;

    vector<int> oth_delta(oth_n);
    vector<int> oth_codes(oth_n);

    for (int i = 0; i < oth_n; i++) {
        in >> oth_delta[i]
           >> oth_codes[i];
    }

    aux.oth_pos = delta_decode(oth_delta);

    aux.oth_ch = decode_other_chars(oth_codes, alphabet);

    return aux;
}

string reconstruct_L3(
    const vector<Match>& matches,
    const vector<Mismatch>& mismatches,
    const vector<int>& ref_encoded
) {
    string result;

    int m_idx = 0;
    int mm_idx = 0;

    while (m_idx < matches.size() ||
           mm_idx < mismatches.size()) {

        if (mm_idx < mismatches.size() &&
            (m_idx >= matches.size() ||
             mismatches[mm_idx].start <
             matches[m_idx].start_target)) {

            result += mismatches[mm_idx].seq;
            mm_idx++;
        }
        else {

            int ref_pos = matches[m_idx].pos;
            int len = matches[m_idx].len;

            for (int i = 0; i < len; i++) {
                result += decode_base(
                    ref_encoded[ref_pos + i]
                );
            }

            m_idx++;
        }
    }

    return result;
}

string reconstruct_L2(
    const string& L3,
    const vector<int>& oth_pos,
    const vector<char>& oth_ch
) {
    string result = L3;

    for (int i = 0; i < (int)oth_pos.size(); i++) {
        result.insert(result.begin() + oth_pos[i], oth_ch[i]);
    }

    return result;
}

string reconstruct_L1(
    const string& L2,
    const vector<int>& N_pos,
    const vector<int>& N_len
) {
    string result = L2;

    for (int i = 0; i < (int)N_pos.size(); i++) {
        result.insert(N_pos[i], string(N_len[i], 'N'));
    }

    return result;
}

string reconstruct_original(
    const string& L1,
    const vector<int>& low_pos,
    const vector<int>& low_len
) {
    string result = L1;

    for (int i = 0; i < (int)low_pos.size(); i++) {
        for (int j = 0; j < low_len[i]; j++) {
            result[low_pos[i] + j] = tolower(result[low_pos[i] + j]);
        }
    }

    return result;
}

void write_fasta(
    const string& filename,
    const string& id,
    const string& sequence,
    const vector<size_t>& seq_len
) {
    ofstream out(filename);

    out << ">" << id << "\n";

    int index = 0;

    for (size_t len : seq_len) {
        out << sequence.substr(index, len) << "\n";
        index += len;
    }
}


void compression_and_decompression_with_output(const std::string& target_file, const std::string& ref_file) {
    try {
        
        FastaData fasta = read_fasta(target_file);

        string genome = join_sequences(fasta.sequences);

        cout << "ID: " << fasta.id << endl;
        cout << "Original genome length: " << genome.size() << endl;

        cout << "\nPREPROCESSING\n" << endl;

        PreprocessedData clean_genome = preprocess_sequence(genome);

        cout << "L1 length: " << clean_genome.L1.size() << endl;
        cout << "L2 length: " << clean_genome.L2.size() << endl;
        cout << "L3 length: " << clean_genome.L3.size() << endl;

        cout << "Lowercase intervals: " << clean_genome.low_pos.size() << endl;

        cout << "N intervals: " << clean_genome.N_pos.size() << endl;

        cout << "Other chars: " << clean_genome.oth_pos.size() << endl;

        cout << "\nREFERENTNI GENOM\n" << endl;

        FastaData fasta_ref = read_fasta(ref_file);

        string genome_ref = join_sequences(fasta_ref.sequences);

        PreprocessedData clean_ref = preprocess_sequence(genome_ref);

        vector<int> ref_encoded = to_binary(clean_ref.L3);

        vector<int> target_encoded = to_binary(clean_genome.L3);

        cout << "Reference encoded length: " << ref_encoded.size() << endl;

        cout << "\nKOMPRESIJA\n" << endl;

        int k = 20;
        int s =  1 << 20;

        vector<Match> matches;
        vector<Mismatch> mismatches;

        greedy_matching(ref_encoded, target_encoded, k, s, matches, mismatches);

        cout << "Matches: " << matches.size() << endl;

        cout << "Mismatches: " << mismatches.size() << endl;

        vector<int> match_pos;

        for (const auto& m : matches) {
            match_pos.push_back(m.pos);
        }

        vector<int> match_pos_delta = delta_encode(match_pos);

        cout << "Delta encoded match positions: " << match_pos_delta.size() << endl;

        vector<int> seq_len_values;
        vector<int> seq_len_counts;

        rle_encode(fasta.seq_len, seq_len_values, seq_len_counts);

        cout << "RLE blocks: " << seq_len_values.size() << endl;

        vector<char> oth_alphabet = build_other_alphabet(clean_genome.oth_ch);

        vector<int> oth_codes = encode_other_chars(clean_genome.oth_ch, oth_alphabet);

        cout << "Other alphabet size: " << oth_alphabet.size() << endl;

        cout << "\n7ZIP KOMPRESIJA\n" << endl;

        write_all_outputs("matches.txt", "mismatches.txt", "auxiliary.txt", matches, mismatches, fasta, clean_genome);

        int zip_status = compress_with_7zip("compressed_output.7z", "matches.txt", "mismatches.txt", "auxiliary.txt");

        if (zip_status != 0) {
            throw runtime_error("7zip compression failed");
        }

        cout << "7zip archive created." << endl;

        cout << "\nDEKOMPRESIJA\n" << endl;

        int extract_status = extract_with_7zip("compressed_output.7z","extracted");

        if (extract_status != 0) {
            throw runtime_error("7zip extraction failed");
        }

        cout << "Archive extracted." << endl;

        vector<Match> decoded_matches = read_matches("extracted/matches.txt");

        cout << "Matches loaded: " << decoded_matches.size() << endl;

        vector<Mismatch> decoded_mismatches = read_mismatches("extracted/mismatches.txt");

        cout << "Mismatches loaded: " << decoded_mismatches.size() << endl;

        AuxData aux = read_auxiliary("extracted/auxiliary.txt");

        cout << "Auxiliary loaded." << endl;

        string reconstructed_L3 = reconstruct_L3(decoded_matches, decoded_mismatches, ref_encoded);

        cout << "L3 reconstructed." << endl;

        string reconstructed_L2 = reconstruct_L2(reconstructed_L3, aux.oth_pos, aux.oth_ch);

        cout << "L2 reconstructed." << endl;

        string reconstructed_L1 = reconstruct_L1(reconstructed_L2, aux.N_pos, aux.N_len);

        cout << "L1 reconstructed." << endl;
        
        string reconstructed_genome = reconstruct_original(reconstructed_L1, aux.low_pos, aux.low_len);

        cout << "Original genome reconstructed." << endl;

        cout << "\nVALIDACIJA\n" << endl;

        if (reconstructed_genome == genome) {
            cout << "USPJEH: Dekompresija je tocna!" << endl;
        }
        else {
            cout << "GRESKA: Sekvence nisu jednake!" << endl;

            int min_len = min(reconstructed_genome.size(), genome.size());

            for (int i = 0; i < min_len; i++) {

                if (reconstructed_genome[i] != genome[i]) {

                    cout << "Prva razlika na poziciji " << i << endl;

                    cout << "Original: " << genome[i] << endl;

                    cout << "Reconstructed: " << reconstructed_genome[i]<< endl;

                    break;
                }
            }
        }
    
        cout << "\nZAPIS REKONSTRUIRANOG FASTA" << endl;

        write_fasta("reconstructed.fna", aux.id, reconstructed_genome, aux.seq_len);

        cout << "Datoteka reconstructed.fna zapisana." << endl;

    }
    catch (const exception& e) {
        cout << "ERROR: " << e.what() << endl;
    }
}

void run_compression_only(const std::string& target_file, const std::string& ref_file){
    FastaData fasta = read_fasta(target_file);
    string genome = join_sequences(fasta.sequences);

    PreprocessedData clean_genome = preprocess_sequence(genome);

    FastaData fasta_ref = read_fasta(ref_file);
    string genome_ref = join_sequences(fasta_ref.sequences);

    PreprocessedData clean_ref = preprocess_sequence(genome_ref);

    vector<int> ref_encoded = to_binary(clean_ref.L3);
    vector<int> target_encoded = to_binary(clean_genome.L3);

    int k = 20;
    int s =  1 << 30;

    vector<Match> matches;
    vector<Mismatch> mismatches;

    greedy_matching(ref_encoded, target_encoded, k, s, matches, mismatches);

    write_all_outputs("matches.txt", "mismatches.txt", "auxiliary.txt", matches, mismatches, fasta, clean_genome);

    int zip_status = compress_with_7zip("compressed_output.7z", "matches.txt", "mismatches.txt", "auxiliary.txt");

    if (zip_status != 0)
        throw runtime_error("7zip compression failed");
}

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cout << "Usage:\n";
        std::cout << "  ./HiRGC output <target.fna> <ref.fna>\n";
        std::cout << "  ./HiRGC totaltime <target.fna> <ref.fna>\n";
        return 1;
    }

    std::string mode = argv[1];
    std::string target_file = argv[2];
    std::string ref_file = argv[3];

    if (mode == "output")
    {
        compression_and_decompression_with_output(target_file, ref_file);
    }
    else if (mode == "totaltime")
    {
        auto start = std::chrono::steady_clock::now();

        run_compression_only(target_file, ref_file);

        auto end = std::chrono::steady_clock::now();

        std::chrono::duration<double> elapsed = end - start;

        std::cout << "Total compression time: "
                  << elapsed.count()
                  << " seconds\n";
    }
    else
    {
        std::cout << "Unknown mode. Use output or benchmark.\n";
        return 1;
    }

    return 0;
}