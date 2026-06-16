#include "runner_mpi.hpp"

int main(int argc, char** argv)
{
  return runner::laderman_mpi(argc, argv, true);
}
