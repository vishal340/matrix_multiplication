# Matrix Multiplication

Fast matrix multiply in C++, with MPI and CUDA variants. Implements Strassen's algorithm and Laderman's 3×3 scheme (23 multiplies instead of 27 per tile). Works on non-square matrices too — they get zero-padded and trimmed on output.

## What's here

**Strassen** — `sequential_strassen`, `strassen_mpi`, `strassen_mpi_no_thread`, `strassen_mpi_2thread`, plus distributed versions (`strassen_mpi_distributed*`) that read binary input. `strassen_cuda` is the CUDA build.

**Laderman** — `laderman_3x3` (timed), `block_3x3` (same thing, no timing), MPI variants (`laderman_mpi*`), and `laderman_cuda` with a real GPU kernel.

**Baseline** — `simple_ikj` for plain O(n³) comparison.

## Build

Needs C++11, CMake 3.10+, and optionally CUDA + MPI depending on what you want to compile.

```bash
mkdir build && cd build
cmake ..
make
```

Binaries land in `build/bin/`.

## Input

Text files: first line is `rows cols`, then the entries row by row (tabs or spaces fine).

```
3 3
1 2 3
4 5 6
7 8 9
```

The distributed Strassen programs want binary matrices instead — use `conv_bin_txt` to convert if you need to inspect them.

## Run

Everything text-based takes three args: `matrix_a matrix_b output`.

```bash
./sequential_strassen matrix_a.txt matrix_b.txt out.txt
./laderman_3x3 matrix_a.txt matrix_b.txt out.txt

mpirun -np 4 ./strassen_mpi matrix_a.txt matrix_b.txt out.txt
mpirun -np 4 ./laderman_mpi matrix_a.txt matrix_b.txt out.txt

./laderman_cuda matrix_a.txt matrix_b.txt out.txt
```

Some builds print timing to stdout before writing the result.

## A few details

Strassen recurses on 2×2 blocks and falls back to naive multiply at 32×32. Laderman tiles everything into 3×3 blocks. `simple_ikj` doesn't pad at all.

MIT license.
