#include "runner_mpi.hpp"

int main(int argc, char** argv)
{
  return runner::strassen_mpi(argc, argv, true);
}
