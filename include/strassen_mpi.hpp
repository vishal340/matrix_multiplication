#ifndef STRASSEN_MPI_HPP
#define STRASSEN_MPI_HPP

#include <memory>

namespace strassen_mpi {

template<typename T>
void add(const T* const A, int jump, const T* const B, int jump1, T* C, int jump2, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      C[i * jump2 + j] = A[i * jump + j] + B[i * jump1 + j];
}

template<typename T>
void atomic_add(const T* const A, int jump, T* B, int jump1, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      B[i * jump1 + j] += A[i * jump + j];
}

template<typename T>
void subtract(const T* const A, int jump, const T* const B, int jump1, T* C, int jump2, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      C[i * jump2 + j] = A[i * jump + j] - B[i * jump1 + j];
}

template<typename T>
void atomic_subtract(const T* const A, int jump, T* B, int jump1, int n1, int n2)
{
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      B[i * jump1 + j] -= A[i * jump + j];
}

template<typename T>
void multiply(const T* const A, int jump, const T* const B, int jump1, T* C, int jump2,
              int n1, int n2, int n3)
{
  for (int j = 0; j < n1; j += 2) {
    for (int i = 0; i < n2; i += 2) {
      for (int k = 0; k < n3; k++) {
        C[j * jump2 + k] += A[i + j * jump] * B[i * jump1 + k];
        C[j * jump2 + k] += A[i + 1 + j * jump] * B[(i + 1) * jump1 + k];
        C[(j + 1) * jump2 + k] += A[i + (j + 1) * jump] * B[i * jump1 + k];
        C[(j + 1) * jump2 + k] += A[(i + 1) + (j + 1) * jump] * B[(i + 1) * jump1 + k];
      }
    }
  }
}

template<typename T>
void C_adjust1(T* C, int jump, int n1, int n2)
{
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      C[i * jump + j] += (C[i * jump + j + n2] - C[(i + n1) * jump + j]);
      C[(i + n1) * jump + j + n2] -= C[i * jump + j];
    }
  }
}

template<typename T>
void C_adjust2(T* C, int jump, int n1, int n2)
{
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      C[(i + n1) * jump + j + n2] += (C[i * jump + j] - C[(i + n1) * jump + j] + C[i * jump + j + n2]);
    }
  }
}

template<typename T>
void strassen(const T* const A, const int jump, const T* const B, const int jump1,
              T* C, const int jump2, int m1, int m2, int m3, int iter)
{
  if (iter == 1)
    multiply(A, jump, B, jump1, C, jump2, m1, m2, m3);
  else {
    iter >>= 1;
    m1 >>= 1;
    m2 >>= 1;
    m3 >>= 1;
    std::unique_ptr<T[]> temp1 = std::make_unique<T[]>(m1 * m2);
    std::unique_ptr<T[]> temp2 = std::make_unique<T[]>(m2 * m3);

    C_adjust1(C, jump2, m1, m3);
    add(&A[0], jump, &A[m2], jump, &temp1[0], m2, m1, m2);
    strassen(&temp1[0], m2, &B[jump1 * m2 + m3], jump1, &C[m3], jump2, m1, m2, m3, iter);
    atomic_subtract(&C[m3], jump2, &C[0], jump2, m1, m3);
    subtract(&B[m3], jump1, &B[m3 + jump1 * m2], jump1, &temp2[0], m3, m2, m3);
    strassen(&A[0], jump, &temp2[0], m3, &C[m3], jump2, m1, m2, m3, iter);
    subtract(&B[jump1 * m2], jump1, &B[0], jump1, &temp2[0], m3, m2, m3);
    strassen(&A[jump * m1 + m2], jump, &temp2[0], m3, &C[jump2 * m1], jump2, m1, m2, m3, iter);
    atomic_add(&C[jump2 * m1], jump2, &C[0], jump2, m1, m3);
    add(&A[jump * m1], jump, &A[jump * m1 + m2], jump, &temp1[0], m2, m1, m2);
    strassen(&temp1[0], m2, &B[0], jump1, &C[jump2 * m1], jump2, m1, m2, m3, iter);
    add(&A[0], jump, &A[jump * m1 + m2], jump, &temp1[0], m2, m1, m2);
    add(&B[0], jump1, &B[jump1 * m2 + m3], jump1, &temp2[0], m3, m2, m3);
    strassen(&temp1[0], m2, &temp2[0], m3, &C[0], jump2, m1, m2, m3, iter);

    C_adjust2(&C[0], jump2, m1, m3);

    subtract(&A[jump * m1], jump, &A[0], jump, &temp1[0], m2, m1, m2);
    add(&B[0], jump1, &B[m3], jump1, &temp2[0], m3, m2, m3);
    strassen(&temp1[0], m2, &temp2[0], m3, &C[jump2 * m1 + m3], jump2, m1, m2, m3, iter);
    subtract(&A[m2], jump, &A[jump * m1 + m2], jump, &temp1[0], m2, m1, m2);
    add(&B[jump1 * m2], jump1, &B[jump1 * m2 + m3], jump1, &temp2[0], m3, m2, m3);
    strassen(&temp1[0], m2, &temp2[0], m3, &C[0], jump2, m1, m2, m3, iter);
  }
}

inline int nearest_ideal(int& n, int& temp, const int temp1, const int threshold)
{
  int t = (temp1 - n % temp1) % temp1;
  int pow = 1;
  n += t;
  int m = (int)(n / temp1);
  while (m > threshold) {
    if (m % 2 == 1) {
      temp += pow;
      m++;
    }
    m >>= 1;
    pow <<= 1;
  }
  temp *= temp1;
  n += temp;
  temp += t;
  if (m % 2 == 1) {
    n += (pow * temp1);
    temp += (pow * temp1);
  }
  return pow;
}

}  // namespace strassen_mpi

#endif
