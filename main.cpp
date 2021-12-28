#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include "IO.cpp"
#include "AffineNeedlemanWunsch.cpp"

typedef std::pair<std::string, std::string> DNA; //<dna_id, dna_string>

void read_dnas(std::istream&, std::vector<DNA>&);

const std::map<char, size_t> default_matrix_map{
{'A', 0},
{'R', 1},
{'N', 2},
{'D', 3},
{'C', 4},
{'Q', 5},
{'E', 6},
{'G', 7},
{'H', 8},
{'I', 9},
{'L', 10},
{'K', 11},
{'M', 12},
{'F', 13},
{'P', 14},
{'S', 15},
{'T', 16},
{'W', 17},
{'Y', 18},
{'V', 19},
{'B', 20},
{'Z', 21},
{'X', 22},
{'*', 23},
};

const char *usage = "Usage: lab3 (-m | --matrix) <input matrix filename> [(-i | --input) <input filename>] [(-o | --output) <output filename>] [--gap-open <gap open>] [--gap-extend <gap extend>]";

int main(int argc, char **argv) {
    std::string input_filename, output_filename, matrix_filename;
    long gap_open = -10;
    long gap_extend = -1;
    
    if (argc % 2 == 0) {
        std::cerr << usage << std::endl;
        return 1;
    }
    
    for (size_t i = 1; i + 1 < argc; i += 2) {
        std::string arg(argv[i]);
        
        if (arg == "-i" || arg == "--input") {
            input_filename = argv[i + 1];
        } else if (arg == "-o" || arg == "--output") {
            output_filename = argv[i + 1];
        } else if (arg == "-gap-open") {
            gap_open = std::stol(argv[i + 1]);
        } else if (arg == "-gap-extend") {
            gap_extend = std::stol(argv[i + 1]);
        } else if (arg == "-m" || arg == "--matrix") {
            matrix_filename = argv[i + 1];
        } else {
            std::cerr << usage << std::endl;
            return 1;
        }
    }
    
    if (matrix_filename.empty()) {
        std::cerr << usage << std::endl;
        return 1;
    }
    
    if (!input_filename.empty()) {
        auto file = std::freopen(input_filename.c_str(), "r", stdin);
        if (!file) {
            std::cerr << "Input file doesn't exist." << std::endl << usage << std::endl;
            return 1;
        }
    }
    
    if (!output_filename.empty()) {
        auto file = std::freopen(output_filename.c_str(), "w", stdout);
        if (!file) {
            std::cerr << "Output file doesn't exist." << std::endl << usage << std::endl;
            return 1;
        }
    }
    
    std::ifstream matrix_ifstream(matrix_filename);
    if (!matrix_ifstream) {
        std::cerr << "Input matrix file doesn't exist." << std::endl << usage << std::endl;
        return 1;
    }
    
    auto matrix = Lobaev::IO::read_matrix<long>(matrix_ifstream);
    
    matrix_ifstream.close();
    
    std::vector<DNA> dnas;
    read_dnas(std::cin, dnas);
    
    if (dnas.size() != 2) {
        std::cerr << "Only 2 dna's in input file are allowed." << std::endl << usage << std::endl;
        return 1;
    }
    
    std::pair<std::pair<std::vector<char>, std::vector<char>>, long> result;
    try {
        result = affine_needleman_wunsch(default_matrix_map, matrix,
                                         std::vector<char>(dnas[0].second.begin(), dnas[0].second.end()),
                                         std::vector<char>(dnas[1].second.begin(), dnas[1].second.end()),
                                         '-',
                                         gap_open,
                                         gap_extend);
    } catch (const std::string &exception) {
        std::cerr << exception << std::endl << usage << std::endl;
        return 1;
    }
    
    const std::string result_string1 = std::string(result.first.first.begin(), result.first.first.end());
    const std::string result_string2 = std::string(result.first.second.begin(), result.first.second.end());
    const long result_score = result.second;
    
    std::cout << result_string1 << std::endl;
    std::cout << result_string2 << std::endl;
    std::cout << "Score: " << result_score << std::endl;
    
    return 0;
}

void read_dnas(std::istream &in, std::vector<DNA> &dnas) {
    std::string cur_line, cur_dna_buf;
    while (std::getline(in, cur_line)) {
        if (cur_line.empty()) {
            continue;
        }
        
        if (cur_line[0] == '>') {
            if (!cur_dna_buf.empty()) {
                dnas.back().second = cur_dna_buf;
                cur_dna_buf.clear();
            }
            
            const size_t id_start_index = cur_line.find('|');
            const size_t id_length = cur_line.substr(id_start_index + 1).find('|');
            const std::string id = cur_line.substr(id_start_index + 1, id_length);
            dnas.emplace_back(id, "");
        } else {
            cur_dna_buf += cur_line;
        }
    }
    
    if (!cur_dna_buf.empty()) {
        dnas.back().second = cur_dna_buf;
    }
}
