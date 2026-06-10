cmake_minimum_required(VERSION 3.10)
project(MatrixMultiplication CUDA CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

# Find required packages
find_package(CUDA REQUIRED)
find_package(MPI REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${MPI_CXX_INCLUDE_DIRS})

# Sequential Strassen implementation
add_executable(sequential_strassen src/sequential_strassen.cpp)

# Simple IKJ multiplication
add_executable(simple_ikj src/simpleIKJ.cpp)

# 3x3 block multiplication
add_executable(block_3x3 src/3_3.cpp)

# MPI Strassen
add_executable(strassen_mpi src/strassen_mpi.cpp)
target_link_libraries(strassen_mpi PUBLIC ${MPI_CXX_LIBRARIES})

# MPI Strassen without threads
add_executable(strassen_mpi_no_thread src/strassen_mpi_without_thread.cpp)
target_link_libraries(strassen_mpi_no_thread PUBLIC ${MPI_CXX_LIBRARIES})

# MPI Strassen with 2 threads
add_executable(strassen_mpi_2thread src/strassen_mpi_2thread.cpp)
target_link_libraries(strassen_mpi_2thread PUBLIC ${MPI_CXX_LIBRARIES} pthread)

# CUDA Strassen
add_executable(strassen_cuda src/strassenRec.cu)
set_target_properties(strassen_cuda PROPERTIES LANGUAGE CUDA)

# Test utilities
add_executable(conv_bin_txt tests/conv_bin_txt.cpp)

# Output directories
set_target_properties(
  sequential_strassen simple_ikj block_3x3 strassen_mpi 
  strassen_mpi_no_thread strassen_mpi_2thread strassen_cuda conv_bin_txt
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
