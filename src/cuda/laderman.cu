#include <iostream>
#include <fstream>
#include <cuda_runtime.h>

#include "laderman.hpp"

using namespace std;

__device__ void multiply_3x3_device(const float a[3][3], const float b[3][3], float c[3][3])
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

__global__ void laderman_multiply_kernel(const float* A, int a_cols, const float* B, int b_cols,
                                         float* C, int c_cols, int n1, int n2, int n3)
{
  int tile_row = blockIdx.y;
  int tile_col = blockIdx.x;
  int i0 = tile_row * 3;
  int j0 = tile_col * 3;
  if (i0 >= n1 || j0 >= n3)
    return;

  float acc[3][3] = {};
  for (int k0 = 0; k0 < n2; k0 += 3) {
    float a[3][3], b[3][3], t[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        a[i][j] = A[(i0 + i) * a_cols + (k0 + j)];
        b[i][j] = B[(k0 + i) * b_cols + (j0 + j)];
      }
    multiply_3x3_device(a, b, t);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        acc[i][j] += t[i][j];
  }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      C[(i0 + i) * c_cols + (j0 + j)] = acc[i][j];
}

void laderman_multiply_gpu(float* A, int a_cols, float* B, int b_cols, float* C, int c_cols,
                           int n1, int n2, int n3)
{
  float *d_A, *d_B, *d_C;
  size_t a_bytes = (size_t)n1 * n2 * sizeof(float);
  size_t b_bytes = (size_t)n2 * n3 * sizeof(float);
  size_t c_bytes = (size_t)n1 * n3 * sizeof(float);

  cudaMalloc(&d_A, a_bytes);
  cudaMalloc(&d_B, b_bytes);
  cudaMalloc(&d_C, c_bytes);
  cudaMemcpy(d_A, A, a_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, B, b_bytes, cudaMemcpyHostToDevice);
  cudaMemset(d_C, 0, c_bytes);

  dim3 grid((n3 + 2) / 3, (n1 + 2) / 3);
  laderman_multiply_kernel<<<grid, 1>>>(d_A, a_cols, d_B, b_cols, d_C, c_cols, n1, n2, n3);
  cudaDeviceSynchronize();

  cudaMemcpy(C, d_C, c_bytes, cudaMemcpyDeviceToHost);
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
}

int main(int argc, char** argv)
{
  ifstream in(argv[1]);
  ifstream in1(argv[2]);
  ofstream out(argv[3]);

  laderman::PaddedMatrices mats = laderman::load_padded_matrices(in, in1);
  in.close();
  in1.close();

  laderman_multiply_gpu(mats.A, mats.n2, mats.B, mats.n3, mats.C, mats.n3,
                        mats.n1, mats.n2, mats.n3);
  laderman::write_result(out, mats);
  laderman::free_matrices(mats);
}
