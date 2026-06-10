# Matrix Multiplication using Strassen's Algorithm

A C++/CUDA program for high-performance matrix multiplication using Strassen's method for both square and non-square matrices.

## Features

- **Sequential Strassen**: CPU-based Strassen algorithm implementation
- **Simple IKJ**: Basic O(n³) matrix multiplication for comparison
- **3x3 Block Multiplication**: Cache-optimized block-based multiplication
- **MPI Strassen**: Distributed memory parallelization using MPI
- **MPI Strassen (2-thread)**: MPI with pthread support for hybrid parallelism
- **CUDA Strassen**: GPU-accelerated implementation

## Project Structure

```
.
├── src/                          # Source implementations
│   ├── sequential_strassen.cpp   # Sequential Strassen algorithm
│   ├── simpleIKJ.cpp             # Simple IKJ multiplication
│   ├── 3_3.cpp                   # 3x3 block multiplication
│   ├── strassen_mpi.cpp          # MPI-based Strassen
│   ├── strassen_mpi_without_thread.cpp
│   ├── strassen_mpi_2thread.cpp  # MPI + threads Strassen
│   └── strassenRec.cu            # CUDA Strassen
├── tests/                        # Test utilities and validation
│   ├── conv_bin_txt.cpp          # Binary to text converter
│   └── test_correct.pp           # Result verification tool
├── include/                      # Header files (future)
├── CMakeLists.txt                # Build configuration
├── .gitignore                    # Git ignore rules
└── README.md                     # This file
```

## Requirements

- C++11 or higher
- CUDA Toolkit (for GPU version)
- OpenMPI or MPICH (for distributed versions)
- CMake 3.10+

## Building

```bash
mkdir build
cd build
cmake ..
make
```

Executables will be generated in `build/bin/`

## Input Format

Matrix files should be formatted as:
```
<rows> <columns>
element_1_1  element_1_2  ...  element_1_c
element_2_1  element_2_2  ...  element_2_c
...
element_r_1  element_r_2  ...  element_r_c
```

Minimum size: 32x32

## Usage Examples

### Sequential Strassen
```bash
./sequential_strassen matrix_a.txt matrix_b.txt output.txt log.txt
```

### MPI Version
```bash
mpirun -np 4 ./strassen_mpi matrix_a.txt matrix_b.txt output.txt log.txt
```

### CUDA Version
```bash
./strassen_cuda matrix_a.txt matrix_b.txt output.txt log.txt
```

## Algorithm Complexity

Strassen's algorithm reduces matrix multiplication complexity from O(n³) to O(n^2.807)

## Verification

Use the test utilities to verify results:
```bash
./conv_bin_txt input_binary.dat input_text.txt
```

## Notes

- The sequential version uses a threshold of 32x32 for switching to standard multiplication
- Non-square matrices are padded to square dimensions with zeros
- Padding values are tracked and removed from output

## License

MIT
