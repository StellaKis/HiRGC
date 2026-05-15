// autori koda: Stella Kiš, Natali Lazarić (cijeli kod je napisan u suradnji)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <sys/resource.h>

using namespace std;

struct FastaData {
    string id; // identifier fasta zapisan nakon znaka '>'
    vector<string> sequences; // vektor koji sadrži sve sekvence iz FASTA datoteke (jedna sekvenca = jedan redak)
    vector<size_t> seq_len; // vektor koji sadrži duljinu svake sekvence
};

struct PreprocessedData {
    string L1;
    string L2;
    string L3;

    vector<int> low_pos; // vektor koji sadrži pozicije malih slova u izvornom zapisu
    vector<int> low_len; // vektor koji sadrži duljine grupa malih slova u izvornom zapisu

    vector<int> N_pos; // vektor koji sadrži pozicije slova 'N' (nepoznatih baza) u izvornom zapisu
    vector<int> N_len; // vektor koji sadrži duljine grupa slova 'N' u izvornom zapisu

    vector<int> oth_pos; // vektor koji sadrži pozicije ostalih znakova u izvornom zapisu
    vector<char> oth_ch; // vektor koji sadrži ostale znakove u izvornom zapisu
};

struct KTupleData {
    vector<uint64_t> values; // vektor koji sadrži vrijednosti k-torki (kodirane kao 64-bitni cijeli brojevi)
    vector<int> positions;  // vektor koji sadrži početne pozicije k-torki u sekvenci (odgovara vrijednostima u 'values')
};

struct HashTable {
    vector<int> h; // vektor koji predstavlja hash tablicu, gdje h[i] je indeks prve k-torke u bucketu
    vector<int> p; // vektor koji predstavlja poveznice između k-torki u istom bucketu, gdje p[i] je indeks sljedeće k-torke u istom bucketu ili -1 ako nema više
};

struct Match {
    int start_target; // početna pozicija podudaranja u ciljnoj sekvenci
    int pos;    // početna pozicija podudaranja u referentnoj sekvenci
    int len;    // duljina podudaranja
};

struct Mismatch {
    int start;  // početna pozicija nepodudaranja u ciljnoj sekvenci
    string seq; // niz znakova koji predstavlja nepodudarni segment u ciljnoj sekvenci
};

struct AuxData {
    string id; // identifier fasta zapisa

    vector<size_t> seq_len; // vektor koji sadrži duljinu svake sekvence (isti kao u FastaData)

    //vektori koji sadrže informacije o pozicijama i duljinama intervala malih slova, slova 'N' i ostalih znakova u izvornom zapisu (isti kao u PreprocessedData)
    vector<int> low_pos; 
    vector<int> low_len; 

    vector<int> N_pos;
    vector<int> N_len;

    vector<int> oth_pos;
    vector<char> oth_ch;
};

// Funkcija za čitanje FASTA datoteke i pohranu podataka u strukturu FastaData
// Ova funkcija popunjava strukturu FastaData s identifikatorom, sekvencama i njihovim duljinama
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


// Funkcija za spajanje vektora sekvenci u jednu cjelovitu sekvencu
string join_sequences(const vector<string>& sequences) {
    string result;

    for (const auto& seq : sequences) {
        result += seq;
    }

    return result;
}

// Funkcija za predobradu sekvence, koja identificira i bilježi intervale malih slova(->L1), slova 'N'(->L2) i ostalih znakova(->L3)
PreprocessedData preprocess_sequence(const string& input) {
    PreprocessedData result;

    bool in_lower = false;
    int start_lower = 0;

    //Identifikazija intervala malih slova
    // bilježenje prvog malog slova u result.low_pos i duljine intervala malih slova u result.low_len
    //i formiranje stringa L1 koji sadrži sve znakove iz inputa pretvorene u velika slova
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
    //Identifikazija intervala slova 'N'
    // bilježenje prvog slova 'N' u result.N_pos i duljine intervala slova 'N' u result.N_len
    //i formiranje stringa L2 koji sadrži sve znakove iz L1 osim slova 'N'

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

    //Identifikazija ostalih znakova (koji nisu A, C, G, T)
    // bilježenje pozicija ostalih znakova u result.oth_pos i samih znakova u result.oth_ch
    //i formiranje stringa L3 koji sadrži samo A, C, G

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

// Funkcija za kodiranje sekvence u binarni format, gdje se A kodira kao 0, C kao 1, G kao 2 i T kao 3
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
// Funkcija za dekodiranje binarnog formata natrag u sekvencu, gdje se 0 dekodira kao A, 1 kao C, 2 kao G i 3 kao T
char decode_base(int x) {
    if (x == 0) return 'A';
    if (x == 1) return 'C';
    if (x == 2) return 'G';
    if (x == 3) return 'T';
    throw runtime_error("Invalid base");
}

// Funkcija za dekodiranje sekvence iz binarnog formata natrag u string, koristeći funkciju decode_base
string decode_sequence(const vector<int>& encoded, int start, int end) {
    string result;

    for (int i = start; i <= end; i++) {
        result += decode_base(encoded[i]);
    }

    return result;
}

// Funkcija za generiranje hash tablice iz k-torki
// Koristi jednostavnu metodu hashiranja (modulo s veličinom tablice) i lančano adresiranje za rješavanje kolizija
HashTable generate_hash_table(const vector<uint64_t>& values, int s) {
    HashTable table;

    int n = (int)values.size();

    table.h.assign(s, -1);  // inicijalno prazno
    table.p.assign(n, -1);  // next pointer

    int mask = s - 1;

    //mapiranje velikog broja (values[i]) u raspon [0, s-1] i izgradnja lančane strukture za k-torke koje se mapiraju u isti bucket
    for (int i = 0; i < n; i++) {
        int hash = (int)(values[i] & mask);

        table.p[i] = table.h[hash]; // p(i) = h(); poveznice između k-torki u istom bucketu
        table.h[hash] = i;  // h() = i; h[j] je indeks prve k-torke u bucketu
    }

    return table;
}

// Funkcija za izvođenje pohlepnog podudaranja između referentne i ciljne sekvence koristeći hash tablicu k-torki
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

    //ako je premalo podataka za formiranje k-torki, cijela ciljna sekvenca se smatra nepodudaranjem (mismatchom)
    if (k <= 0 || nr < k || nt < k) {
        if (nt > 0) {
            Mismatch mm;
            mm.start = 0;
            mm.seq = decode_sequence(target_encoded, 0, nt - 1);
            mismatches.push_back(mm);
        }
        return;
    }

    //generiranje k-torki iz referentne sekvence
    KTupleData ref_tuples;
    int ref_tuple_count = nr - k + 1; //broj svih k-torki u referenci

    //kodiranje prve k-torke referentne sekvence u 64-bitni cijeli broj
    uint64_t ref_value = 0;

    for (int j = k - 1; j >= 0; j--) {
        ref_value <<= 2;
        ref_value += ref_encoded[j];
    }

    ref_tuples.values.push_back(ref_value);
    ref_tuples.positions.push_back(0);

    //kodiranje preostalih k-torki referentne sekvence koristeći klizni prozor i bitovne operacije
    int shift_bits = 2 * (k - 1);

    for (int i = 1; i < ref_tuple_count; i++) {
        ref_value >>= 2; // uklanjanje najstarijeg elementa k-torke 
        ref_value += ((uint64_t)ref_encoded[i + k - 1] << shift_bits); // dodavanje novog elementa k-torke na najviši bit

        ref_tuples.values.push_back(ref_value); // pohrana nove kodirane k-torke
        ref_tuples.positions.push_back(i); // pohrana početne pozicije nove k-torke u referentnoj sekvenci
    }

    //generiranje hash tablice iz k-torki referentne sekvence
    HashTable table = generate_hash_table(ref_tuples.values, s);

    //hash lookup i pohlepno podudaranje
    int mask = s - 1;

    int i = 0; //trenutna pozicija u ciljnoj sekvenci
    int p_star = 0; //pozicija od koje počinje sljedeći mismatch u ciljnoj sekvenci tj. početak zadnjeg neobrađenog dijela (za mismatch)

    // računa se prva k-torka ciljne sekvence (kao broj u bazi 4)
    uint64_t Vt = 0; // hash vrijednost trenutne k-torke ciljne sekvence

    for (int j = k - 1; j >= 0; j--) {
        Vt <<= 2;
        Vt += target_encoded[j];
    }

    int current_hash_pos = 0; //pozicija u ciljnoj sekvenci do koje je već izračunata hash vrijednost (Vt)

    while (i <= nt - k) {
        while (current_hash_pos < i) {
            Vt >>= 2; // uklanjanje najstarijeg elementa k-torke
            Vt += ((uint64_t)target_encoded[current_hash_pos + k] << shift_bits); // dodavanje novog elementa k-torke na najviši bit
            current_hash_pos++; // pomicanje pozicije do koje je izračunat hash
        }

        int bucket = (int)(Vt & mask); // određivanje bucket-a u hash tablici za trenutnu k-torku ciljne sekvence
        int j = table.h[bucket]; // indeks prve k-torke u referentnoj sekvenci koja se mapira u isti bucket

        int p_max = -1; // pozicija u referentnoj sekvenci do koje se podudara trenutna k-torka
        int l_max = 0; // duljina podudaranja između trenutne k-torke ciljne sekvence i odgovarajuće k-torke referentne sekvence

        while (j != -1) { 
            //ako hash vrijednost k-torke iz referentne sekvence ne odgovara hash vrijednosti trenutne k-torke ciljne sekvence, nastavlja se pretraživanje sljedeće k-torke u istom bucketu
            if (ref_tuples.values[j] != Vt) {
                j = table.p[j];
                continue;
            }

            //ako hash vrijednosti odgovaraju, provjerava se stvarno podudaranje k-torki i produžuje se podudaranje koliko god je moguće
            int ref_pos = ref_tuples.positions[j]; // početna pozicija k-torke u referentnoj sekvenci

            int l = k; // početna duljina podudaranja (barem k, jer se hash vrijednosti poklapaju)

            // produžavanje podudaranja dok se znakovi poklapaju i ne dođe do kraja jedne od sekvenci
            while ((i + l < nt) &&
                (ref_pos + l < nr) &&
                target_encoded[i + l] == ref_encoded[ref_pos + l]) {
                l++;
            }

            // ako je pronađeno duže podudaranje od dosadašnjeg maksimalnog, ažurira se maksimalna duljina i pozicija podudaranja u referentnoj sekvenci
            if (l > l_max) {
                p_max = ref_pos;
                l_max = l;
            }

            // nastavlja se pretraživanje sljedeće k-torke u istom bucketu
            j = table.p[j];
        }

        //ako je pronađeno podudaranje (l_max > 0), bilježi se podudaranje i eventualno prethodno nepodudaranje (ako postoji) u odgovarajuće vektore matches i mismatches
        if (l_max > 0) {
            //ako trenutna pozicija u ciljnoj sekvenci (i) je veća od p_star, znači da postoji nepodudaranje između p_star i i-1, koje se bilježi kao mismatch
            if (i > p_star) {
                Mismatch mm;
                mm.start = p_star;
                mm.seq = decode_sequence(target_encoded, p_star, i - 1);
                mismatches.push_back(mm);
            }

            //bilježi se podudaranje u vektor matches s početnom pozicijom u ciljnoj sekvenci, pozicijom u referentnoj sekvenci i duljinom podudaranja
            Match m;
            m.start_target = i;
            m.pos = p_max;
            m.len = l_max;
            matches.push_back(m);

            p_star = i + l_max; // ažurira se p_star na poziciju nakon kraja pronađenog podudaranja, što označava početak sljedećeg neobrađenog dijela ciljne sekvence
        }

        i = i + l_max + 1; // pomiče se trenutna pozicija u ciljnoj sekvenci na kraj pronađenog podudaranja + 1, što znači da se sljedeće pretraživanje započinje odmah nakon kraja pronađenog podudaranja
    }

    // ako nakon završetka glavne petlje postoji neobrađeni dio ciljne sekvence (od p_star do kraja), bilježi se kao mismatch
    if (p_star < nt) {
        Mismatch mm;
        mm.start = p_star;
        mm.seq = decode_sequence(target_encoded, p_star, nt - 1); 
        mismatches.push_back(mm);
    }
}

//delta kodiranje = prva vrijednost ostaje nepromijenjena, a svaka sljedeća zapisuje se kao razlika u odnosu na prethodnu
vector<int> delta_encode(const vector<int>& values) {
    vector<int> result;

    if (values.empty()) return result;

    result.push_back(values[0]); // prva vrijednost ostaje nepromijenjena

    for (int i = 1; i < (int)values.size(); i++) {
        result.push_back(values[i] - values[i - 1]); // svaka sljedeća vrijednost zapisuje se kao razlika u odnosu na prethodnu
    }

    return result;
}

// funkcija za run-length encoding (RLE) služi za kodiranje duljina redaka kako bi se mogla obnoviti izvorna struktura FASTA datoteke
void rle_encode(const vector<size_t>& input, vector<int>& values, vector<int>& counts) {
    values.clear();
    counts.clear();

    if (input.empty()) return;

    int current = (int)input[0]; // trenutna vrijednost koja se kodira (na početku prva duljina reda)
    int count = 1; // brojač ponavljanja trenutne vrijednosti

    for (int i = 1; i < (int)input.size(); i++) {
        if ((int)input[i] == current) { // ako je trenutna vrijednost jednaka idućoj (input[i]), povećava se brojač ponavljanja i provjerava za idući i
            count++;
        } else { // ako se vrijednost ne podudara 
            values.push_back(current); //bilježi se trenutna vrijednost
            counts.push_back(count); //i njen broj ponavljanja
            current = (int)input[i]; //zatim se ažurira trenutna vrijednost 
            count = 1; //i resetira brojač
        }
    }

    values.push_back(current);
    counts.push_back(count);
}

// funkcija za izgradnju abecede ostalih znakova koji nisu A, C, G, T, na temelju vektora ostalih znakova (oth_ch) prikupljenih tijekom predobrade
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

// funkcija za kodiranje ostalih znakova (koji nisu A, C, G, T) u cijele brojeve na temelju izgrađene abecede ostalih znakova
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

// funkcija za pisanje svih izlaza (matches, mismatches i pomoćnih podataka) u odgovarajuće datoteke koje će se zatim komprimirati
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

    { //pohrana delta kodiranih matches u datoteku
        ofstream out(matches_file);
        if (!out.is_open()) throw runtime_error("Cannot open matches file");

        out << matches.size() << "\n";
        for (int i = 0; i < (int)matches.size(); i++) {
            out << match_start_delta[i] << " " << match_pos_delta[i] << " " << match_len[i] << "\n";
        }
    }

    { //pohrana mismatches u datoteku
        ofstream out(mismatches_file);
        if (!out.is_open()) throw runtime_error("Cannot open mismatches file");

        out << mismatches.size() << "\n";
        for (const auto& mm : mismatches) {
            out << mm.start << " " << mm.seq << "\n";
        }
    }

    {   //pohrana FASTA header-a (identifikatora), (RLE enkodirane) duljine redaka, (delta enkodiranih) informacija o malim slovima,  N intervala i ostale znakove u auxiliary datoteku
        ofstream out(aux_file);
        if (!out.is_open()) throw runtime_error("Cannot open auxiliary file");

        out << "ID\n" << fasta.id << "\n";
        // pohrana RLE enkodiranih duljina redaka
        out << "SEQ_LEN_RLE\n" << seq_len_values.size() << "\n";
        for (int i = 0; i < (int)seq_len_values.size(); i++) {
            out << seq_len_values[i] << " " << seq_len_counts[i] << "\n";
        }
        // pohrana delta enkodiranih informacija o malim slovima
        out << "LOWERCASE\n" << prep.low_pos.size() << "\n";
        for (int i = 0; i < (int)prep.low_pos.size(); i++) {
            out << low_pos_delta[i] << " " << prep.low_len[i] << "\n";
        }
        // pohrana delta enkodiranih informacija o N intervalima
        out << "N_INTERVALS\n" << prep.N_pos.size() << "\n";
        for (int i = 0; i < (int)prep.N_pos.size(); i++) {
            out << N_pos_delta[i] << " " << prep.N_len[i] << "\n";
        }
        // pohrana abecede ostalih znakova 
        out << "OTHER_ALPHABET\n" << oth_alphabet.size() << "\n";
        for (char c : oth_alphabet) {
            out << c << " ";
        }
        out << "\n";
        // pohrana delta enkodiranih pozicija ostalih znakova i njihovih kodova u abecedi
        out << "OTHER_CHARS\n";
        out << prep.oth_pos.size() << "\n";
        for (int i = 0; i < (int)prep.oth_pos.size(); i++) {
            out << oth_pos_delta[i] << " " << oth_codes[i] << "\n";
        }
    }
}

// komprimiranje izlaznih datoteka koristeći 7zip s PPMd algoritmom i određenom razinom kompresije = 5
int compress_with_7zip(const string& archive_name, const string& matches_file, const string& mismatches_file, const string& aux_file) {
    string command = "7z a -t7z -m0=PPMd -mx=5 \"" + archive_name + "\" \"" + matches_file + "\" \"" + mismatches_file + "\" \"" + aux_file + "\"";

    return system(command.c_str());
}

//dekompresija
// funkcija za dekompresiju arhive koristeći 7zip, koja će se koristiti u funkciji HiRGC dekompresije 
int extract_with_7zip(const string& archive_name, const string& output_dir) {
    string command = "7z x \"" + archive_name + "\" -o\"" + output_dir + "\" -y";

    return system(command.c_str());
}

// funkcija za dekodiranje delta kodiranih vrijednosti natrag u originalne vrijednosti, koristeći činjenicu da je prva vrijednost nepromijenjena, a svaka sljedeća se dobiva zbrajanjem prethodne dekodirane vrijednosti i trenutne delta vrijednosti
vector<int> delta_decode(const vector<int>& delta) {
    vector<int> result;

    if (delta.empty()) return result;

    result.push_back(delta[0]);

    for (int i = 1; i < (int)delta.size(); i++) {
        result.push_back(result[i - 1] + delta[i]); //svaka sljedeća vrijednost se dobiva zbrajanjem prethodne dekodirane vrijednosti i trenutne delta vrijednosti
    }

    return result;
}

// funkcija za dekodiranje RLE kodiranih vrijednosti natrag u originalne vrijednosti, dobivanje duljine redaka iz RLE zapisa
vector<size_t> rle_decode(const vector<int>& values, const vector<int>& counts) {
    vector<size_t> result;

    // svaka vrijednost se ponavlja onoliko puta koliko je navedeno u counts, i dodaje se u rezultat
    for (int i = 0; i < (int)values.size(); i++) {
        for (int j = 0; j < counts[i]; j++) {
            result.push_back(values[i]);
        }
    }

    return result;
}

// funkcija za dekodiranje ostalih znakova (koji nisu A, C, G, T) natrag u originalne znakove koristeći kodove i abecedu ostalih znakova
vector<char> decode_other_chars(const vector<int>& codes, const vector<char>& alphabet) {
    vector<char> result;

    for (int code : codes) {
        result.push_back(alphabet[code]);
    }

    return result;
}

// funkcija za čitanje matches iz datoteka koje su prethodno dekomprimirane, kako bi se mogla rekonstruirati originalna sekvenca
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
    // dekodiranje delta kodiranih početnih pozicija i pozicija podudaranja natrag u originalne vrijednosti
    vector<int> start = delta_decode(start_delta);
    vector<int> pos = delta_decode(pos_delta);

    vector<Match> matches(n);
    // popunjavanje vektora matches s dekodiranim vrijednostima početnih pozicija, pozicija podudaranja i duljina podudaranja
    for (int i = 0; i < n; i++) {
        matches[i].start_target = start[i];
        matches[i].pos = pos[i];
        matches[i].len = len[i];
    }

    return matches;
}

// funkcija za čitanje mismatches iz datoteka koje su prethodno dekomprimirane, kako bi se mogla rekonstruirati originalna sekvenca
vector<Mismatch> read_mismatches(const string& filename) {
    ifstream in(filename);

    if (!in.is_open()) {
        throw runtime_error("Cannot open mismatches");
    }

    int n;
    in >> n;

    vector<Mismatch> mismatches(n);
    // popunjavanje vektora mismatches s početnim pozicijama i nepodudarnim segmentima pročitanim iz datoteke
    for (int i = 0; i < n; i++) {
        in >> mismatches[i].start
           >> mismatches[i].seq;
    }

    return mismatches;
}

// funkcija za čitanje pomoćnih podataka iz auxiliary datoteke koja je prethodno dekomprimirana, kako bi se mogla rekonstruirati originalna sekvenca
AuxData read_auxiliary(const string& filename) {

    ifstream in(filename);

    if (!in.is_open()) {
        throw runtime_error("Cannot open auxiliary file");
    }

    AuxData aux;

    string token;

    // pohrana FASTA header-a (identifikatora) u pomoćnu strukturu AuxData

    in >> token; // ID
    in.ignore();

    getline(in, aux.id);

    in >> token;

    int rle_n;
    in >> rle_n;

    vector<int> values(rle_n);
    vector<int> counts(rle_n);
    // pohrana RLE kodiranih duljina redaka iz auxilary datoteke u vektore values i counts
    for (int i = 0; i < rle_n; i++) {
        in >> values[i] >> counts[i];
    }
    // dekodiranje RLE kodiranih vrijednosti natrag u originalne vrijednosti tj duljine redaka i pohrana u pomoćnu strukturu AuxData
    aux.seq_len = rle_decode(values, counts);

    //LOWERCASE

    in >> token;

    int low_n;
    in >> low_n;

    vector<int> low_delta(low_n);

    aux.low_len.resize(low_n);
    // pohrana delta kodiranih početnih pozicija intervala malih slova i njihovih duljina iz auxiliary datoteke u vektore low_delta i aux.low_len
    for (int i = 0; i < low_n; i++) {
        in >> low_delta[i]
           >> aux.low_len[i];
    }
    // dekodiranje delta kodiranih početnih pozicija intervala malih slova natrag u originalne vrijednosti i pohrana u pomoćnu strukturu AuxData
    aux.low_pos = delta_decode(low_delta);

    //N_INTERVALS

    in >> token;

    int N_n;
    in >> N_n;

    vector<int> N_delta(N_n);

    aux.N_len.resize(N_n);
    // pohrana delta kodiranih početnih pozicija intervala slova 'N' i njihovih duljina iz auxiliary datoteke u vektore N_delta i aux.N_len
    for (int i = 0; i < N_n; i++) {
        in >> N_delta[i]
           >> aux.N_len[i];
    }
    // dekodiranje delta kodiranih početnih pozicija intervala slova 'N' natrag u originalne vrijednosti i pohrana u pomoćnu strukturu AuxData
    aux.N_pos = delta_decode(N_delta);

    //OTHER_ALPHABET

    in >> token;

    int alphabet_size;
    in >> alphabet_size;

    vector<char> alphabet(alphabet_size);
    // pohrana abecede ostalih znakova koji nisu A, C, G, T iz auxiliary datoteke u vektor alphabet
    for (int i = 0; i < alphabet_size; i++) {
        in >> alphabet[i];
    }

    //OTHER_CHARS

    in >> token;

    int oth_n;
    in >> oth_n;

    vector<int> oth_delta(oth_n);
    vector<int> oth_codes(oth_n);
    // pohrana delta kodiranih pozicija ostalih znakova i njihovih kodova iz auxiliary datoteke u vektore oth_delta i oth_codes
    for (int i = 0; i < oth_n; i++) {
        in >> oth_delta[i]
           >> oth_codes[i];
    }
    // dekodiranje delta kodiranih pozicija ostalih znakova natrag u originalne vrijednosti i pohrana u pomoćnu strukturu AuxData
    aux.oth_pos = delta_decode(oth_delta);
    // dekodiranje kodova ostalih znakova natrag u originalne znakove koristeći izgrađenu abecedu ostalih znakova i pohrana u pomoćnu strukturu AuxData
    aux.oth_ch = decode_other_chars(oth_codes, alphabet);

    return aux;
}

// funkcija za rekonstrukciju stringa L3 koristeći informacije o podudaranjima i nepodudaranjima te kodiranu referentnu sekvencu
string reconstruct_L3( const vector<Match>& matches, const vector<Mismatch>& mismatches, const vector<int>& ref_encoded) {
    string result;

    int m_idx = 0;
    int mm_idx = 0;
    //

    while (m_idx < matches.size() ||
           mm_idx < mismatches.size()) {
        //postoji li još nepodudaranja i pojavljuje li se trenutno nepodudaranje prije sljedećeg podudaranja ?
        if (mm_idx < mismatches.size() &&
            (m_idx >= matches.size() ||
             mismatches[mm_idx].start <
             matches[m_idx].start_target)) {
            // ako da, dodaje se nepodudaranje u rezultat i pomiče se indeks nepodudaranja
            result += mismatches[mm_idx].seq;
            mm_idx++;
        }
        else {
            //ako ne, dodaje se podudaranje u rezultat i pomiče se indeks podudaranja
            int ref_pos = matches[m_idx].pos;
            int len = matches[m_idx].len;
            // dekodiranje podudaranja iz referentne sekvence koristeći kodiranu referentnu sekvencu (ref_encoded) i funkciju decode_base, te dodavanje dekodiranog podudaranja u rezultat
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

// funkcija za rekonstrukciju stringa L2 dodavanjem ostalih znakova na njihove originalne pozicije u stringu L3
string reconstruct_L2(const string& L3, const vector<int>& oth_pos, const vector<char>& oth_ch) {
    string result = L3;

    for (int i = 0; i < (int)oth_pos.size(); i++) {
        result.insert(result.begin() + oth_pos[i], oth_ch[i]);
    }

    return result;
}

// funkcija za rekonstrukciju stringa L1 dodavanjem intervala slova 'N' na njihove originalne pozicije u stringu L2
string reconstruct_L1(const string& L2, const vector<int>& N_pos, const vector<int>& N_len) {
    string result = L2;

    for (int i = 0; i < (int)N_pos.size(); i++) {
        result.insert(N_pos[i], string(N_len[i], 'N'));
    }

    return result;
}

// funkcija za rekonstrukciju originalne sekvence dodavanjem intervala malih slova na njihove originalne pozicije u stringu L1
string reconstruct_original(const string& L1, const vector<int>& low_pos, const vector<int>& low_len) {
    string result = L1;

    for (int i = 0; i < (int)low_pos.size(); i++) {
        for (int j = 0; j < low_len[i]; j++) {
            result[low_pos[i] + j] = tolower(result[low_pos[i] + j]);
        }
    }

    return result;
}

// funkcija za pisanje rekonstruirane sekvence u FASTA datoteku, koristeći informacije o originalnim duljinama redaka kako bi se sačuvala struktura FASTA datoteke
void write_fasta( const string& filename, const string& id, const string& sequence, const vector<size_t>& seq_len) {
    ofstream out(filename);
    // zapisivanje FASTA header-a (identifikatora) u datoteku
    out << ">" << id << "\n";

    int index = 0;
    // pisanje sekvence u datoteku u redcima čije duljine su zadane u vektoru seq_len, kako bi se sačuvala struktura FASTA datoteke
    for (size_t len : seq_len) {
        out << sequence.substr(index, len) << "\n";
        index += len;
    }
}

// funkcija za mjerenje utroška memorije
long get_total_peak_memory_kb() {
    struct rusage self_usage;
    struct rusage children_usage;

    getrusage(RUSAGE_SELF, &self_usage);
    getrusage(RUSAGE_CHILDREN, &children_usage);

    return self_usage.ru_maxrss + children_usage.ru_maxrss;
}

// funkcija koja demonstrira cijeli proces kompresije i dekompresije, uključujući čitanje sekvenci, predobradu, pohlepno podudaranje, kodiranje, pisanje izlaza, kompresiju s 7zipom, dekompresiju i rekonstrukciju originalne sekvence
void compression_and_decompression_with_output(const std::string& target_file, const std::string& ref_file) {
    try {
        
        // Čitanje sekvenci iz FASTA datoteke
        FastaData fasta = read_fasta(target_file);
        // Spajanje svih sekvenci u jednu cjelovitu sekvencu
        string genome = join_sequences(fasta.sequences);

        cout << "ID: " << fasta.id << endl;
        cout << "Original genome length: " << genome.size() << endl;

        cout << "\nPREPROCESSING\n" << endl;

        // Predobrada sekvence: formiranje stringova L1, L2 i L3 (identifikacija intervala malih slova, slova 'N' i ostalih znakova)
        PreprocessedData clean_genome = preprocess_sequence(genome);

        cout << "L1 length: " << clean_genome.L1.size() << endl;
        cout << "L2 length: " << clean_genome.L2.size() << endl;
        cout << "L3 length: " << clean_genome.L3.size() << endl;

        cout << "Lowercase intervals: " << clean_genome.low_pos.size() << endl;

        cout << "N intervals: " << clean_genome.N_pos.size() << endl;

        cout << "Other chars: " << clean_genome.oth_pos.size() << endl;

        cout << "\nREFERENTNI GENOM\n" << endl;
        // Čitanje referentne sekvence iz FASTA datoteke i predobrada na isti način kao i ciljne sekvence
        FastaData fasta_ref = read_fasta(ref_file);

        string genome_ref = join_sequences(fasta_ref.sequences);

        PreprocessedData clean_ref = preprocess_sequence(genome_ref);

        vector<int> ref_encoded = to_binary(clean_ref.L3);

        vector<int> target_encoded = to_binary(clean_genome.L3);

        cout << "Reference encoded length: " << ref_encoded.size() << endl;

        cout << "\nKOMPRESIJA\n" << endl;
        // Postavljanje parametara k i s za pohlepno podudaranje
        int k = 20;
        int s =  1 << 20;

        vector<Match> matches;
        vector<Mismatch> mismatches;
        // Izvođenje pohlepnog podudaranja između referentne i ciljne sekvence, bilježenje podudaranja i nepodudaranja
        greedy_matching(ref_encoded, target_encoded, k, s, matches, mismatches);

        cout << "Matches: " << matches.size() << endl;

        cout << "Mismatches: " << mismatches.size() << endl;

        vector<int> match_pos;

        for (const auto& m : matches) {
            match_pos.push_back(m.pos);
        }
        // Delta kodiranje pozicija podudaranja
        vector<int> match_pos_delta = delta_encode(match_pos);

        cout << "Delta encoded match positions: " << match_pos_delta.size() << endl;

        vector<int> seq_len_values;
        vector<int> seq_len_counts;

        // RLE kodiranje duljina sekvenci
        rle_encode(fasta.seq_len, seq_len_values, seq_len_counts);

        cout << "RLE blocks: " << seq_len_values.size() << endl;
        // Izgradnja abecede ostalih znakova i kodiranje tih znakova u brojeve
        vector<char> oth_alphabet = build_other_alphabet(clean_genome.oth_ch);

        vector<int> oth_codes = encode_other_chars(clean_genome.oth_ch, oth_alphabet);

        cout << "Other alphabet size: " << oth_alphabet.size() << endl;

        cout << "\n7ZIP KOMPRESIJA\n" << endl;
        // Zapis svih izlaza (matches, mismatches, auxiliary) u tekstualne datoteke i kompresija tih datoteka u 7zip arhivu
        write_all_outputs("matches.txt", "mismatches.txt", "auxiliary.txt", matches, mismatches, fasta, clean_genome);

        int zip_status = compress_with_7zip("compressed_output.7z", "matches.txt", "mismatches.txt", "auxiliary.txt");

        if (zip_status != 0) {
            throw runtime_error("7zip compression failed");
        }

        cout << "7zip archive created." << endl;

        cout << "\nDEKOMPRESIJA\n" << endl;
        // Ekstrakcija datoteka iz 7zip arhive, čitanje podataka iz tih datoteka i rekonstrukcija originalne sekvence kroz obrnut proces od kompresije
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

// funkcija koja pokreće samo proces dekompresije,
void decompression (const std::string& target_file, const std::string& ref_file){

        FastaData fasta_ref = read_fasta(ref_file);

        string genome_ref = join_sequences(fasta_ref.sequences);

        PreprocessedData clean_ref = preprocess_sequence(genome_ref);

        vector<int> ref_encoded = to_binary(clean_ref.L3);

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

        FastaData fasta = read_fasta(target_file);

        string genome = join_sequences(fasta.sequences);

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

// funkcija koja pokreće samo proces kompresije, bez dekompresije i validacije, kako bi se moglo izmjeriti ukupno vrijeme kompresije
void run_compression_only(const std::string& target_file, const std::string& ref_file){
    // Čitanje sekvenci iz FASTA datoteke i spajanje svih sekvenci u jednu cjelovitu sekvencu
    FastaData fasta = read_fasta(target_file);
    string genome = join_sequences(fasta.sequences);
    // Predobrada sekvence: formiranje stringova L1, L2 i L3 (identifikacija intervala malih slova, slova 'N' i ostalih znakova)
    PreprocessedData clean_genome = preprocess_sequence(genome);
    // Čitanje referentne sekvence iz FASTA datoteke i predobrada na isti način kao i ciljne sekvence
    FastaData fasta_ref = read_fasta(ref_file);
    string genome_ref = join_sequences(fasta_ref.sequences);

    PreprocessedData clean_ref = preprocess_sequence(genome_ref);
    // Kodiranje sekvenci L3 referentne i ciljne sekvence u binarni format, gdje se A kodira kao 0, C kao 1, G kao 2 i T kao 3
    vector<int> ref_encoded = to_binary(clean_ref.L3);
    vector<int> target_encoded = to_binary(clean_genome.L3);

    int k = 20;
    int s =  1 << 30;

    vector<Match> matches;
    vector<Mismatch> mismatches;
    // Izvođenje pohlepnog podudaranja između referentne i ciljne sekvence, bilježenje podudaranja i nepodudaranja
    greedy_matching(ref_encoded, target_encoded, k, s, matches, mismatches);
    // Zapis svih izlaza (matches, mismatches, auxiliary) u tekstualne datoteke i kompresija tih datoteka u 7zip arhivu
    write_all_outputs("matches.txt", "mismatches.txt", "auxiliary.txt", matches, mismatches, fasta, clean_genome);
    int zip_status = compress_with_7zip("compressed_output.7z", "matches.txt", "mismatches.txt", "auxiliary.txt");

    if (zip_status != 0)
        throw runtime_error("7zip compression failed");
}

// glavna funkcija koja omogućuje pokretanje programa u dva načina: s izlazom (kompresija i dekompresija s validacijom) ili samo s mjerenjem ukupnog vremena kompresije
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
        
        long memory_kb = get_total_peak_memory_kb();

        std::cout << "Peak total memory usage: "
                << memory_kb / 1024.0 << " MB\n";
    }
    else if (mode == "decompression")
    {
        decompression(target_file, ref_file);

    }
    else
    {
        std::cout << "Unknown mode. Use output or totaltime.\n";
        return 1;
    }

    return 0;
}