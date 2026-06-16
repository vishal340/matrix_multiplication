cmake_minimum_required(VERSION 3.10)

# Generate compile_commands.json for clangd and other language servers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile commands" FORCE)

project(MatrixMultiplication CUDA CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

find_package(CUDA REQUIRED)
find_package(MPI REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${MPI_CXX_INCLUDE_DIRS})

# --- Core libraries ---
add_library(strassen_core src/core/strassen.cpp)
target_include_directories(strassen_core PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(laderman_core src/core/laderman.cpp)
target_include_directories(laderman_core PUBLIC ${CMAKE_SOURCE_DIR}/include)

# --- Baseline ---
add_executable(simple_ikj src/baseline/simpleIKJ.cpp)

# --- Strassen ---
add_executable(sequential_strassen src/strassen/sequential.cpp)
target_link_libraries(sequential_strassen PRIVATE strassen_core)

add_executable(strassen_mpi src/mpi/strassen.cpp)
target_link_libraries(strassen_mpi PUBLIC ${MPI_CXX_LIBRARIES} strassen_core)

add_executable(strassen_mpi_no_thread src/mpi/strassen_timed.cpp)
target_link_libraries(strassen_mpi_no_thread PUBLIC ${MPI_CXX_LIBRARIES} strassen_core)

add_executable(strassen_mpi_2thread src/mpi/strassen_2thread.cpp)
target_link_libraries(strassen_mpi_2thread PUBLIC ${MPI_CXX_LIBRARIES} strassen_core pthread)

add_executable(strassen_mpi_distributed src/mpi/distributed/strassen.cpp)
target_link_libraries(strassen_mpi_distributed PUBLIC ${MPI_CXX_LIBRARIES})

add_executable(strassen_mpi_distributed_no_thread src/mpi/distributed/strassen_timed.cpp)
target_link_libraries(strassen_mpi_distributed_no_thread PUBLIC ${MPI_CXX_LIBRARIES})

add_executable(strassen_mpi_distributed_2thread src/mpi/distributed/strassen_2thread.cpp)
target_link_libraries(strassen_mpi_distributed_2thread PUBLIC ${MPI_CXX_LIBRARIES} pthread)

add_executable(strassen_cuda src/cuda/strassen.cu)
set_target_properties(strassen_cuda PROPERTIES LANGUAGE CUDA)

# --- Laderman ---
add_executable(laderman_3x3 src/laderman/sequential.cpp)
target_link_libraries(laderman_3x3 PRIVATE laderman_core)

add_executable(block_3x3 src/laderman/block.cpp)
target_link_libraries(block_3x3 PRIVATE laderman_core)

add_executable(laderman_mpi src/mpi/laderman.cpp)
target_link_libraries(laderman_mpi PUBLIC ${MPI_CXX_LIBRARIES} laderman_core)

add_executable(laderman_mpi_no_thread src/mpi/laderman_timed.cpp)
target_link_libraries(laderman_mpi_no_thread PUBLIC ${MPI_CXX_LIBRARIES} laderman_core)

add_executable(laderman_mpi_2thread src/mpi/laderman_2thread.cpp)
target_link_libraries(laderman_mpi_2thread PUBLIC ${MPI_CXX_LIBRARIES} laderman_core pthread)

add_executable(laderman_cuda src/cuda/laderman.cu)
set_target_properties(laderman_cuda PROPERTIES LANGUAGE CUDA)

# --- Utilities ---
add_executable(conv_bin_txt tests/conv_bin_txt.cpp)

set(ALL_TARGETS
  simple_ikj
  sequential_strassen strassen_mpi strassen_mpi_no_thread strassen_mpi_2thread
  strassen_mpi_distributed strassen_mpi_distributed_no_thread strassen_mpi_distributed_2thread
  strassen_cuda
  laderman_3x3 block_3x3 laderman_mpi laderman_mpi_no_thread laderman_mpi_2thread
  laderman_cuda
  conv_bin_txt
)
set_target_properties(${ALL_TARGETS} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(CMAKE_EXPORT_COMPILE_COMMANDS AND NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  set(_compile_commands_link "${CMAKE_SOURCE_DIR}/compile_commands.json")
  if(EXISTS "${_compile_commands_link}")
    file(REMOVE "${_compile_commands_link}")
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E create_symlink
            "${CMAKE_BINARY_DIR}/compile_commands.json"
            "${_compile_commands_link}"
    RESULT_VARIABLE _compile_commands_symlink_status
    ERROR_QUIET
  )
  if(NOT _compile_commands_symlink_status EQUAL 0)
    message(STATUS "compile_commands.json will be at ${CMAKE_BINARY_DIR}/compile_commands.json")
  endif()
endif()
