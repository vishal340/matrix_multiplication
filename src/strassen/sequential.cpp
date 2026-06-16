#include "runner.hpp"

int main(int argc, char** argv)
{
  return runner::strassen_sequential(argc, argv, true);
}
