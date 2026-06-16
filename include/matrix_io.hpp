#ifndef MATRIX_IO_HPP
#define MATRIX_IO_HPP

#include <algorithm>
#include <fstream>

#include "strassen.hpp"

namespace matrix_io {

constexpr int block3 = 3;

inline int pad_to_block(int& n, int& temp, int block)
{
  temp = (block - n % block) % block;
  n += temp;
  return block;
}

struct PaddedMatrices {
  float* A = nullptr;
  float* B = nullptr;
  float* C = nullptr;
  int n1 = 0;
  int n2 = 0;
  int n3 = 0;
  int temp1 = 0;
  int temp2 = 0;
  int temp3 = 0;
  int iter = 0;
};

inline void write_unpadded(std::ofstream& out, const float* C, int c_cols, const PaddedMatrices& m)
{
  out << m.n1 << "\t" << m.n3 << "\n";
  for (int i = 0; i < m.n1 - m.temp1; i++) {
    for (int j = 0; j < m.n3 - m.temp3; j++)
      out << C[i * c_cols + j] << "\t";
    out << "\n";
  }
}

inline void free_matrices(PaddedMatrices& m)
{
  delete[] m.C;
  delete[] m.A;
  delete[] m.B;
}

inline PaddedMatrices load_strassen(std::ifstream& in, std::ifstream& in1)
{
  PaddedMatrices m{};
  in >> m.n1 >> m.n2;
  in1 >> m.n2 >> m.n3;

  m.iter = strassen::nearest_ideal(m.n1, m.temp1);
  m.iter = std::min(m.iter, strassen::nearest_ideal(m.n2, m.temp2));
  m.iter = std::min(m.iter, strassen::nearest_ideal(m.n3, m.temp3));

  m.A = new float[m.n1 * m.n2];
  m.B = new float[m.n2 * m.n3];
  m.C = new float[m.n1 * m.n3]();

  for (int i = 0; i < m.n1 - m.temp1; i++) {
    for (int j = 0; j < m.n2 - m.temp2; j++)
      in >> m.A[i * m.n2 + j];
    for (int j = m.n2 - m.temp2; j < m.n2; j++)
      m.A[i * m.n2 + j] = 0;
  }
  for (int i = (m.n1 - m.temp1) * m.n2; i < m.n1 * m.n2; i++)
    m.A[i] = 0;

  for (int i = 0; i < m.n2 - m.temp2; i++) {
    for (int j = 0; j < m.n3 - m.temp3; j++)
      in1 >> m.B[i * m.n3 + j];
    for (int j = m.n3 - m.temp3; j < m.n3; j++)
      m.B[i * m.n3 + j] = 0;
  }
  for (int i = (m.n2 - m.temp2) * m.n3; i < m.n2 * m.n3; i++)
    m.B[i] = 0;

  return m;
}

inline PaddedMatrices load_laderman(std::ifstream& in, std::ifstream& in1)
{
  PaddedMatrices m{};
  in >> m.n1 >> m.n2;
  in1 >> m.n2 >> m.n3;

  pad_to_block(m.n1, m.temp1, block3);
  pad_to_block(m.n2, m.temp2, block3);
  pad_to_block(m.n3, m.temp3, block3);

  m.A = new float[m.n1 * m.n2]();
  m.B = new float[m.n2 * m.n3]();
  m.C = new float[m.n1 * m.n3]();

  for (int i = 0; i < m.n1 - m.temp1; i++)
    for (int j = 0; j < m.n2 - m.temp2; j++)
      in >> m.A[i * m.n2 + j];

  for (int i = 0; i < m.n2 - m.temp2; i++)
    for (int j = 0; j < m.n3 - m.temp3; j++)
      in1 >> m.B[i * m.n3 + j];

  return m;
}

inline void write_result(std::ofstream& out, const PaddedMatrices& m)
{
  write_unpadded(out, m.C, m.n3, m);
}

inline PaddedMatrices load_padded_matrices(std::ifstream& in, std::ifstream& in1)
{
  return load_strassen(in, in1);
}

}  // namespace matrix_io

#endif
