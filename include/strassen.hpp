#ifndef STRASSEN_HPP
#define STRASSEN_HPP

#include <cstring>

namespace strassen {

constexpr int threshold = 32;

inline void add(float* A, int jump, float* B, int jump1, float* C, int jump2, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      C[i * jump2 + j] = A[i * jump + j] + B[i * jump1 + j];
}

inline void atomic_add(float* A, int jump, float* B, int jump1, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      B[i * jump1 + j] += A[i * jump + j];
}

inline void subtract(float* A, int jump, float* B, int jump1, float* C, int jump2, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      C[i * jump2 + j] = A[i * jump + j] - B[i * jump1 + j];
}

inline void atomic_subtract(float* A, int jump, float* B, int jump1, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      B[i * jump1 + j] -= A[i * jump + j];
}

inline void multiply(float* A, int jump, float* B, int jump1, float* C, int jump2,
                     int n1, int n2, int n3)
{
  for (int j = 0; j < n1; j += 2)
    for (int i = 0; i < n2; i += 2)
      for (int k = 0; k < n3; k++) {
        C[j * jump2 + k] += A[i + j * jump] * B[i * jump1 + k];
        C[j * jump2 + k] += A[i + 1 + j * jump] * B[(i + 1) * jump1 + k];
        C[(j + 1) * jump2 + k] += A[i + (j + 1) * jump] * B[i * jump1 + k];
        C[(j + 1) * jump2 + k] += A[(i + 1) + (j + 1) * jump] * B[(i + 1) * jump1 + k];
      }
}

inline void nonzero_C(float* C, int jump2, int m1, int m3)
{
  atomic_add(C + m3 / 2, jump2, C, jump2, m1 / 2, m3 / 2);
  atomic_add(C + jump2 * m1 / 2, jump2, C + jump2 * m1 / 2 + m3 / 2, jump2, m1 / 2, m3 / 2);
  atomic_subtract(C, jump2, C + jump2 * m1 / 2 + m3 / 2, jump2, m1 / 2, m3 / 2);
}

void strassen(float* A, int jump, float* B, int jump1, float* C, int jump2,
              int n1, int n2, int n3, int iter, bool flag);

inline int nearest_ideal(int& n, int& temp)
{
  temp = 0;
  int pow = 1;
  int m = n;
  while (m > threshold) {
    if (m % 2 == 1) {
      temp += pow;
      m++;
    }
    m /= 2;
    pow *= 2;
  }
  n += temp;
  return pow;
}

}  // namespace strassen

#endif
