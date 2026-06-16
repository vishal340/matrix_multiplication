#ifndef RUNNER_HPP
#define RUNNER_HPP

#include <fstream>
#include <iostream>
#include <sys/time.h>

#include "laderman.hpp"
#include "matrix_io.hpp"
#include "strassen.hpp"

namespace runner {

inline void print_elapsed(const timespec& start, const timespec& end)
{
  std::cout << "Time taken: "
            << (end.tv_nsec - start.tv_nsec) + static_cast<double>(end.tv_sec - start.tv_sec) * 1e9;
}

inline int strassen_sequential(int argc, char** argv, bool timed)
{
  if (argc < 4)
    return 1;

  std::ifstream in(argv[1]);
  std::ifstream in1(argv[2]);
  std::ofstream out(argv[3]);
  matrix_io::PaddedMatrices mats = matrix_io::load_strassen(in, in1);
  in.close();
  in1.close();

  timespec start{}, end{};
  if (timed)
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

  strassen::strassen(mats.A, mats.n2, mats.B, mats.n3, mats.C, mats.n3,
                     mats.n1, mats.n2, mats.n3, mats.iter, 1);

  if (timed) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    print_elapsed(start, end);
  }

  matrix_io::write_result(out, mats);
  matrix_io::free_matrices(mats);
  return 0;
}

inline int laderman_sequential(int argc, char** argv, bool timed)
{
  if (argc < 4)
    return 1;

  std::ifstream in(argv[1]);
  std::ifstream in1(argv[2]);
  std::ofstream out(argv[3]);
  laderman::PaddedMatrices mats = laderman::load_padded_matrices(in, in1);
  in.close();
  in1.close();

  timespec start{}, end{};
  if (timed)
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

  laderman::multiply(mats.A, mats.n2, mats.B, mats.n3, mats.C, mats.n3,
                     mats.n1, mats.n2, mats.n3);

  if (timed) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    print_elapsed(start, end);
  }

  laderman::write_result(out, mats);
  laderman::free_matrices(mats);
  return 0;
}

}  // namespace runner

#endif
