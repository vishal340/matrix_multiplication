#ifndef RUNNER_MPI_HPP
#define RUNNER_MPI_HPP

#include <mpi.h>
#include <pthread.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/time.h>

#include "laderman.hpp"
#include "matrix_io.hpp"
#include "runner.hpp"
#include "strassen.hpp"

namespace runner {

inline int strassen_mpi(int argc, char** argv, bool timed)
{
  if (argc < 4)
    return 1;

  MPI_Init(&argc, &argv);
  std::ifstream in(argv[1]);
  std::ifstream in1(argv[2]);
  std::ofstream out(argv[3]);
  matrix_io::PaddedMatrices mats = matrix_io::load_strassen(in, in1);
  in.close();
  in1.close();

  timespec start{}, end{};
  if (timed) {
    MPI_Barrier(MPI_COMM_WORLD);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  }

  strassen::strassen(mats.A, mats.n2, mats.B, mats.n3, mats.C, mats.n3,
                     mats.n1, mats.n2, mats.n3, mats.iter, 1);

  if (timed) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    MPI_Barrier(MPI_COMM_WORLD);
    print_elapsed(start, end);
  }

  matrix_io::write_result(out, mats);
  matrix_io::free_matrices(mats);
  MPI_Finalize();
  return 0;
}

inline int laderman_mpi(int argc, char** argv, bool timed)
{
  if (argc < 4)
    return 1;

  MPI_Init(&argc, &argv);
  std::ifstream in(argv[1]);
  std::ifstream in1(argv[2]);
  std::ofstream out(argv[3]);
  laderman::PaddedMatrices mats = laderman::load_padded_matrices(in, in1);
  in.close();
  in1.close();

  timespec start{}, end{};
  if (timed) {
    MPI_Barrier(MPI_COMM_WORLD);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  }

  laderman::multiply(mats.A, mats.n2, mats.B, mats.n3, mats.C, mats.n3,
                     mats.n1, mats.n2, mats.n3);

  if (timed) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    MPI_Barrier(MPI_COMM_WORLD);
    print_elapsed(start, end);
  }

  laderman::write_result(out, mats);
  laderman::free_matrices(mats);
  MPI_Finalize();
  return 0;
}

struct LadermanThreadParam {
  float* A;
  int a_cols;
  float* B;
  int b_cols;
  float* C;
  int c_cols;
  int n2;
  int n3;
  int row_start;
  int row_end;
};

inline void* laderman_thread_entry(void* arg)
{
  LadermanThreadParam* param = static_cast<LadermanThreadParam*>(arg);
  laderman::multiply_row_range(param->A, param->a_cols, param->B, param->b_cols,
                               param->C, param->c_cols, param->n2, param->n3,
                               param->row_start, param->row_end);
  free(param);
  return NULL;
}

inline int laderman_mpi_2thread(int argc, char** argv)
{
  if (argc < 4)
    return 1;

  MPI_Init(&argc, &argv);
  std::ifstream in(argv[1]);
  std::ifstream in1(argv[2]);
  std::ofstream out(argv[3]);
  laderman::PaddedMatrices mats = laderman::load_padded_matrices(in, in1);
  in.close();
  in1.close();

  timespec start{}, end{};
  MPI_Barrier(MPI_COMM_WORLD);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

  pthread_t thread1, thread2;
  LadermanThreadParam* param1 =
      static_cast<LadermanThreadParam*>(malloc(sizeof(LadermanThreadParam)));
  LadermanThreadParam* param2 =
      static_cast<LadermanThreadParam*>(malloc(sizeof(LadermanThreadParam)));
  int mid = (mats.n1 / laderman::block_size / 2) * laderman::block_size;

  param1->A = mats.A;
  param1->a_cols = mats.n2;
  param1->B = mats.B;
  param1->b_cols = mats.n3;
  param1->C = mats.C;
  param1->c_cols = mats.n3;
  param1->n2 = mats.n2;
  param1->n3 = mats.n3;
  param1->row_start = 0;
  param1->row_end = mid;

  param2->A = mats.A;
  param2->a_cols = mats.n2;
  param2->B = mats.B;
  param2->b_cols = mats.n3;
  param2->C = mats.C;
  param2->c_cols = mats.n3;
  param2->n2 = mats.n2;
  param2->n3 = mats.n3;
  param2->row_start = mid;
  param2->row_end = mats.n1;

  pthread_create(&thread1, NULL, laderman_thread_entry, param1);
  pthread_create(&thread2, NULL, laderman_thread_entry, param2);
  pthread_join(thread1, NULL);
  pthread_join(thread2, NULL);

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
  MPI_Barrier(MPI_COMM_WORLD);
  print_elapsed(start, end);

  laderman::write_result(out, mats);
  laderman::free_matrices(mats);
  MPI_Finalize();
  return 0;
}

}  // namespace runner

#endif
