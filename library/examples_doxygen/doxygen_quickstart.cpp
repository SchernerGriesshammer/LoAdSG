/**
 * \page quickstart_page Quick Start
 *
# Installation

1. Clone repository
   ```sh
   git clone
   ```
2. create build folder
   ```sh
   cd expdesg/library
   mkdir build
   cd build
   ```
3. Run CMake
   ```sh
   cmake -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="./main.cpp" ..
   ```
   For options see the Options section below.

5. make the programm
   ```sh
   make
   ```
   Use ``make -j8`` to use 8 threads.

6. Run program
   ```sh
   ./sgrun
   ```

7. Run parallel
   ```sh
   mpirun -np numberofnodes ./sgrun
   ```
### Options
 * ``-DCMAKE_BUILD_TYPE={Release,Debug}``
 * ``-DMAIN_FILE="./main2.cpp"`` with target main file in reference to the ``library`` directory
 * ``-DMPI_ON={true,false}`` for compiling with mpi
 * ``-DOMP_ON={true,false}`` for compiling with omp
 *

  Note that using MPI and/or OMP requires setting up a valid environment.  See \ref usage_example_page   "Usage Examples".


### Example calls
   ```sh
   cmake -DMPI_ON=false -DOMP_ON=false -DCUDA_ON=false -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="./main.cpp" .. && make -j4 && ./sgrun
   ```

# Change main_file
### Using CMake argument
Call cmake with ``-DMAIN_FILE="./main2.cpp"``:
   ```sh
   cmake -DMPI_ON=false -DOMP_ON=false -DCMAKE_BUILD_TYPE=Debug -DMAIN_FILE="./main2.cpp" ..
   ```
### Or by overwriting default in CMakeList.txt
   In
   ```sh
   ../expdesg/library/CmakeLists.txt
   ```
   replace main_file.cpp with desired main file:
   ```sh
   add_executable(${PROJECT_NAME} main_file.cpp)
   ```
 *
 *
 *
 */

