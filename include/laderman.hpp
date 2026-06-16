#ifndef LADERMAN_HPP
#define LADERMAN_HPP

#include "matrix_io.hpp"

namespace laderman {

constexpr int block_size = matrix_io::block3;

void multiply_3x3(const float a[3][3], const float b[3][3], float c[3][3]);

void multiply_3x3_accumulate(const float* A, int a_stride, const float* B, int b_stride,
                             float* C, int c_stride);

void multiply(float* A, int a_cols, float* B, int b_cols, float* C, int c_cols,
              int n1, int n2, int n3);

void multiply_row_range(float* A, int a_cols, float* B, int b_cols, float* C, int c_cols,
                        int n1, int n2, int n3, int row_start, int row_end);

using PaddedMatrices = matrix_io::PaddedMatrices;

inline PaddedMatrices load_padded_matrices(std::ifstream& in, std::ifstream& in1)
{
  return matrix_io::load_laderman(in, in1);
}

inline void write_result(std::ofstream& out, const PaddedMatrices& m)
{
  matrix_io::write_result(out, m);
}

inline void free_matrices(PaddedMatrices& m)
{
  matrix_io::free_matrices(m);
}

}  // namespace laderman

#endif
