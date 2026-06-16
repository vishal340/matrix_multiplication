#include "strassen.hpp"

namespace strassen {

void strassen(float* A, int jump, float* B, int jump1, float* C, int jump2,
              int n1, int n2, int n3, int iter, bool flag)
{
  if (iter == 1)
    multiply(A, jump, B, jump1, C, jump2, n1, n2, n3);
  else {
    iter /= 2;
    int m1 = n1 / 2, m2 = n2 / 2, m3 = n3 / 2;
    float* temp1 = new float[m1 * m2];
    float* temp2 = new float[m2 * m3];

    add(A, jump, A + jump * m1 + m2, jump, temp1, m2, m1, m2);
    add(B, jump1, B + jump1 * m2 + m3, jump1, temp2, m3, m2, m3);
    if (iter > 1 && !flag)
      nonzero_C(C, jump2, m1, m3);
    strassen(temp1, m2, temp2, m3, C, jump2, m1, m2, m3, iter, flag);
    if (flag) {
      for (int i = 0; i < m1; i++)
        for (int j = 0; j < m3; j++)
          C[jump2 * m1 + m3 + i * jump2 + j] = C[i * jump2 + j];
    } else
      atomic_add(C, jump2, C + jump2 * m1 + m3, jump2, m1, m3);

    add(A + jump * m1, jump, A + jump * m1 + m2, jump, temp1, m2, m1, m2);
    if (iter > 1 && !flag)
      nonzero_C(C + jump2 * m1, jump2, m1, m3);
    strassen(temp1, m2, B, jump1, C + jump2 * m1, jump2, m1, m2, m3, iter, flag);

    add(A, jump, A + m2, jump, temp1, m2, m1, m2);
    if (iter > 1 && !flag)
      nonzero_C(C + m3, jump2, m1, m3);
    strassen(temp1, m2, B + jump1 * m2 + m3, jump1, C + m3, jump2, m1, m2, m3, iter, flag);
    atomic_subtract(C + m3, jump2, C, jump2, m1, m3);

    subtract(A + jump * m1, jump, A, jump, temp1, m2, m1, m2);
    add(B, jump1, B + m3, jump1, temp2, m3, m2, m3);

    if (iter > 1)
      nonzero_C(C + jump2 * m1 + m3, jump2, m1, m3);
    strassen(temp1, m2, temp2, m3, C + jump2 * m1 + m3, jump2, m1, m2, m3, iter, flag);

    atomic_subtract(C + jump2 * m1, jump2, C + jump2 * m1 + m3, jump2, m1, m3);

    subtract(A + m2, jump, A + jump * m1 + m2, jump, temp1, m2, m1, m2);
    add(B + jump1 * m2, jump1, B + jump1 * m2 + m3, jump1, temp2, m3, m2, m3);

    if (iter > 1)
      nonzero_C(C, jump2, m1, m3);
    strassen(temp1, m2, temp2, m3, C, jump2, m1, m2, m3, iter, 0);

    delete[] temp1;
    temp1 = new float[m1 * m3]();

    subtract(B + m3, jump1, B + m3 + jump1 * m2, jump1, temp2, m3, m2, m3);
    strassen(A, jump, temp2, m3, temp1, m3, m1, m2, m3, iter, 1);
    atomic_add(temp1, m3, C + m3, jump2, m1, m3);
    atomic_add(temp1, m3, C + jump2 * m1 + m3, jump2, m1, m3);

    memset(temp1, 0, sizeof(float) * m1 * m3);
    subtract(B + jump1 * m2, jump1, B, jump1, temp2, m3, m2, m3);
    strassen(A + jump * m1 + m2, jump, temp2, m3, temp1, m3, m1, m2, m3, iter, 1);
    atomic_add(temp1, m3, C + jump2 * m1, jump2, m1, m3);
    atomic_add(temp1, m3, C, jump2, m1, m3);
    delete[] temp1;
    delete[] temp2;
  }
}

}  // namespace strassen
