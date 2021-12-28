#include "Matrix.cpp"

namespace Lobaev::Math {
    
    template <class T, class L>
    std::pair<std::pair<std::vector<T>, std::vector<T>>, L> affine_needleman_wunsch(const std::map<T, size_t> &matrix_map,
                                                                                    const Matrix<L> &matrix,
                                                                                    const std::vector<T> &seq1,
                                                                                    const std::vector<T> &seq2,
                                                                                    const T seq_gap,
                                                                                    const L gap_open,
                                                                                    const L gap_extend
    ) {
        if (matrix.rows_count() != matrix_map.size() || matrix.columns_count() != matrix_map.size()) {
            throw "affine_needleman_wunsch: matrix size should be " +
                  std::to_string(matrix_map.size()) + "*" + std::to_string(matrix_map.size());
        }
        
        const L infinity = 2 * gap_open + (seq1.size() + seq2.size()) * gap_extend + 1;
    
        Matrix<L> matrix_m(seq1.size() + 1, seq2.size() + 1);
        for (size_t row_index = 1; row_index < matrix_m.rows_count(); row_index++) {
            matrix_m(row_index, 0) = infinity;
        }
        for (size_t column_index = 1; column_index < matrix_m.columns_count(); column_index++) {
            matrix_m(0, column_index) = infinity;
        }
        
        Matrix<L> matrix_i(seq1.size() + 1, seq2.size() + 1);
        for (size_t row_index = 0; row_index < matrix_i.rows_count(); row_index++) {
            matrix_i(row_index, 0) = infinity;
        }
        for (size_t column_index = 1; column_index < matrix_i.columns_count(); column_index++) {
            matrix_i(0, column_index) = gap_open + (column_index - 1) * gap_extend;
        }
        
        Matrix<L> matrix_d(seq1.size() + 1, seq2.size() + 1);
        for (size_t row_index = 1; row_index < matrix_d.rows_count(); row_index++) {
            matrix_d(row_index, 0) = gap_open + (row_index - 1) * gap_extend;
        }
        for (size_t column_index = 0; column_index < matrix_d.columns_count(); column_index++) {
            matrix_d(0, column_index) = infinity;
        }
        
        for (size_t i = 1; i <= seq1.size(); i++) {
            for (size_t j = 1; j <= seq2.size(); j++) {
                const size_t map_index_i = matrix_map.at(seq1[i - 1]);
                const size_t map_index_j = matrix_map.at(seq2[j - 1]);
                
                matrix_m(i, j) = std::max(matrix_m(i - 1, j - 1),
                                          std::max(matrix_i(i - 1, j - 1), matrix_d(i - 1, j - 1))) +
                                          matrix(map_index_i, map_index_j);
                
                matrix_i(i, j) = std::max(matrix_i(i, j - 1) + gap_extend,
                                          std::max(matrix_m(i, j - 1), matrix_d(i, j - 1)) + gap_open);
                
                matrix_d(i, j) = std::max(matrix_d(i - 1, j) + gap_extend,
                                          std::max(matrix_m(i - 1, j), matrix_i(i - 1, j)) + gap_open);
            }
        }
        
        std::pair<std::vector<T>, std::vector<T>> result_sequences;
        size_t i = seq1.size(), j = seq2.size();
        while (i > 0 && j > 0) {
            const L cur_max_value = std::max(matrix_m(i, j),
                                             std::max(matrix_i(i, j),
                                                      matrix_d(i, j)));
            if (cur_max_value == matrix_d(i, j)) {
                i--;
                result_sequences.first.push_back(seq1[i]);
                result_sequences.second.push_back(seq_gap);
            } else if (cur_max_value == matrix_i(i, j)) {
                j--;
                result_sequences.first.push_back(seq_gap);
                result_sequences.second.push_back(seq2[j]);
            } else {
                i--;
                j--;
                result_sequences.first.push_back(seq1[i]);
                result_sequences.second.push_back(seq2[j]);
            }
        }
        std::reverse(result_sequences.first.begin(), result_sequences.first.end());
        std::reverse(result_sequences.second.begin(), result_sequences.second.end());
        
        return std::make_pair(result_sequences, std::max(matrix_m(seq1.size(), seq2.size()),
                                                         std::max(matrix_i(seq1.size(), seq2.size()),
                                                                  matrix_d(seq1.size(), seq2.size()))));
    }
    
}
