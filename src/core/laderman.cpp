#include "laderman.hpp"

namespace laderman {

// Julian Laderman (1976): 3x3 matrix product in 23 multiplications.
void multiply_3x3(const float a[3][3], const float b[3][3], float c[3][3])
{
  float m[24];

  m[1] = (a[0][0] + a[0][1] + a[0][2] - a[1][0] - a[1][1] - a[2][1] - a[2][2]) * b[1][1];
  m[2] = (a[0][0] - a[1][0]) * (-b[0][1] + b[1][1]);
  m[3] = a[1][1] * (-b[0][0] + b[0][1] + b[1][0] - b[1][1] - b[1][2] - b[2][0] + b[2][2]);
  m[4] = (-a[0][0] + a[1][0] + a[1][1]) * (b[0][0] - b[0][1] + b[1][1]);
  m[5] = (a[1][0] + a[1][1]) * (-b[0][0] + b[0][1]);
  m[6] = a[0][0] * b[0][0];
  m[7] = (-a[0][0] + a[2][0] + a[2][1]) * (b[0][0] - b[0][2] + b[1][2]);
  m[8] = (-a[0][0] + a[2][0]) * (b[0][2] - b[1][2]);
  m[9] = (a[2][0] + a[2][1]) * (-b[0][0] + b[0][2]);
  m[10] = (a[0][0] + a[0][1] + a[0][2] - a[1][1] - a[1][2] - a[2][0] - a[2][1]) * b[1][2];
  m[11] = a[2][1] * (-b[0][0] + b[0][2] + b[1][0] - b[1][1] - b[1][2] - b[2][0] + b[2][1]);
  m[12] = (-a[0][2] + a[2][1] + a[2][2]) * (b[1][1] + b[2][0] - b[2][1]);
  m[13] = (a[0][2] - a[2][2]) * (b[1][1] - b[2][1]);
  m[14] = a[0][2] * b[2][0];
  m[15] = (a[2][1] + a[2][2]) * (-b[2][0] + b[2][1]);
  m[16] = (-a[0][2] + a[1][1] + a[1][2]) * (b[1][2] + b[2][0] - b[2][2]);
  m[17] = (a[0][2] - a[1][2]) * (b[1][2] - b[2][2]);
  m[18] = (a[1][1] + a[1][2]) * (-b[2][0] + b[2][2]);
  m[19] = a[0][1] * b[1][0];
  m[20] = a[1][2] * b[2][1];
  m[21] = a[1][0] * b[0][2];
  m[22] = a[2][0] * b[0][1];
  m[23] = a[2][2] * b[2][2];

  c[0][0] = m[6] + m[14] + m[19];
  c[0][1] = m[1] + m[4] + m[5] + m[6] + m[12] + m[14] + m[15];
  c[0][2] = m[6] + m[7] + m[9] + m[10] + m[14] + m[16] + m[18];
  c[1][0] = m[2] + m[3] + m[4] + m[6] + m[14] + m[16] + m[17];
  c[1][1] = m[2] + m[4] + m[5] + m[6] + m[20];
  c[1][2] = m[14] + m[16] + m[17] + m[18] + m[21];
  c[2][0] = m[6] + m[7] + m[8] + m[11] + m[12] + m[13] + m[14];
  c[2][1] = m[12] + m[13] + m[14] + m[15] + m[22];
  c[2][2] = m[6] + m[7] + m[8] + m[9] + m[23];
}

void multiply_3x3_accumulate(const float* A, int a_stride, const float* B, int b_stride,
                             float* C, int c_stride)
{
  float a[3][3], b[3][3], t[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      a[i][j] = A[i * a_stride + j];
      b[i][j] = B[i * b_stride + j];
    }
  multiply_3x3(a, b, t);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      C[i * c_stride + j] += t[i][j];
}

void multiply(float* A, int a_cols, float* B, int b_cols, float* C, int c_cols,
              int n1, int n2, int n3)
{
  multiply_row_range(A, a_cols, B, b_cols, C, c_cols, n1, n2, n3, 0, n1);
}

void multiply_row_range(float* A, int a_cols, float* B, int b_cols, float* C, int c_cols,
                        int n1, int n2, int n3, int row_start, int row_end)
{
  for (int i0 = row_start; i0 < row_end; i0 += block_size)
    for (int j0 = 0; j0 < n3; j0 += block_size)
      for (int k0 = 0; k0 < n2; k0 += block_size)
        multiply_3x3_accumulate(A + i0 * a_cols + k0, a_cols,
                                  B + k0 * b_cols + j0, b_cols,
                                  C + i0 * c_cols + j0, c_cols);
}

}  // namespace laderman
