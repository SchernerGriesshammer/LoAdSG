# Installation

1. Clone repo
   ```sh
   git clone 
   ```
2. create build forlder
   ```sh
   cd library
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="./main2.cpp" ..
   make
   ```
3. Run CMake
   ```sh
   cmake -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="./main2.cpp" ..
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
## Options
 * ``-DCMAKE_BUILD_TYPE={Release,Debug}``
 * ``-DMAIN_FILE="./main2.cpp"`` with target main file in reference to the ``library`` directory
 * ``-DMPI_ON={true,false}`` for compiling with mpi
 * ``-DOMP_ON={true,false}`` for compiling with omp
 * ``-DCUDA_ON={true,false}`` for compiling with CUDA
  
  Note that using MPI, OMP and CUDA together requires setting up a valid environment.
* MPI+CUDA: it is possible to set the number of MPI tasks that should be allocated per GPU. 

## Example calls
   ```sh
   cmake -DMPI_ON=false -DOMP_ON=false -DCUDA_ON=false -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="./poisson_neumann_test.cpp" .. && make -j4 && ./sgrun 
   ```
# Change main_file
## Using CMake argument
Call cmake with ``-DMAIN_FILE="./main2.cpp"``:
   ```sh
   cmake -DMPI_ON=false -DOMP_ON=false -DCMAKE_BUILD_TYPE=Debug -DMAIN_FILE="./main2.cpp" ..
   ```
## Or by overwriting default in CMakeList.txt
   In
   ```sh
   ../LoAdSG/library/CmakeLists.txt
   ```
   replace main_file.cpp with desired main file:
   ```sh
   add_executable(${PROJECT_NAME} main_file.cpp)
   ```

# Generate Performance Graph
### GPerfTools
Note that GPerfTools are not supported by WSL2.

Install 
```
sudo apt-get install google-perftools libgoogle-perftools-dev
```
Call the programm:
```
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so CPUPROFILE=./text.prof ./sgrun
```
Get the output as pdf:
```
google-pprof -pdf ./sgrun text.prof  > ./text.pdf
```

### Perf

For method runtimes use:
```
perf record ./sgrun
perf report
```

For Hardware resources use:

```
 perf stat -d ./sgrun
```

### LIKWID
```
 likwid-perfctr -C 1 -g CACHES ./sgrun
 ```
# RRZE
## MEGGIE
* Get available Modules: `` module avail ``
* Load Modules: `` module load <name> ``
* Get Interactive Job: `` salloc -N1 ``
* Start Program: go to ``scripts`` folder and call: ``./RunJob.sh <MPI> <FILE> <NODES> <CORES> gcc`` example ``./RunJob.sh 2 main2_test.cpp 2 40 gcc``
  * On Meggie there are 20 cores per Node. Allocating more Cores than the Number of Nodes *20 will result in an crash
  * For help call  ``./RunJob.sh``

## TinyGPU
```sh
salloc -N1 --gres=gpu:rtx2080ti:{numGPUs} --cpus-per-task={OMP_Threads per MPI thread} --ntasks={numMPI_Threads} --partition=work

module load cmake/3.18.4 && module load cuda/11.6.1 && module load intelmpi/2019.8

cmake -DMPI_ON=true -DOMP_ON=true -DCUDA_ON=true -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="poisson_neumann_test.cpp" .. && make -j8

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun ./sgrun
```
## Fritz
 cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DMPI_ON=false -DOMP_ON=false -DCUDA_ON=false -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="poisson_neumann_test.cpp" .. && make -j8 && ./sgrun  

# CUDA
cmake -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_CXX_COMPILER=g++-8 -DMPI_ON=false -DOMP_ON=false -DCUDA_ON=true -DCMAKE_BUILD_TYPE=Release -DMAIN_FILE="poisson_neumann_test.cpp" .. && make -j8 && ./sgrun

# RRZE setup

### SSH Access per key
* `cd ~/.ssh`
* ``ssh-keygen`` here I named the file ``rrze_key``
* ``ssh-copy-id -i rrze_key.pub <accname>@cshpc.rrze.fau.de``
* add to the ``config`` file inside the .ssh folder following lines:
  ```ssh
  Host cshpc
      HostName cshpc.rrze.fau.de
      IdentityFile ~/.ssh/rrze_key
      User <accname>
  Host meggie
      HostName meggie.rrze.fau.de
      IdentityFile ~/.ssh/rrze_key
      User <accname>
      ProxyJump cshpc
  ```
* now you can connect to meggie with ``ssh meggie``
* alternatively connect to cshpc with ``ssh cshpc`` 
  
### Setup gitlab use
on cshpc or meggie or any other rrze system do:
* `cd ~/.ssh`
* ``ssh-keygen`` here I named the file ``gitlab_key``
* now ```vim gitlab_key.pub`` and copy the contents
* now go to the gitlab website "https://gitlab.cs.fau.de/-/profile/keys" and under SSH key paste the contents and press "add key"
* now in the terminal do ```vim ~/.ssh/config`` and add
    ```ssh
  Host gitlab.cs.fau.de
      HostName gitlab.cs.fau.de
      IdentityFile ~/.ssh/gitlab_key
      User root
  ``` 
