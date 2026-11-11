#!/usr/bin/env bash
set -euo pipefail

library_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${library_dir}/build"

rm -rf "${build_dir}"

cmake -S "${library_dir}" -B "${build_dir}" \
    -DMAIN_FILE=mainfiles/eigen_solver_benchmark.cpp \
    -DMPI_ON=true \
    -DOPENMP_ON=true \
    -DCMAKE_BUILD_TYPE=Release

cmake --build "${build_dir}" -j"$(nproc)"

cd "${build_dir}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
    mpirun -np "${MPI_PROCESSES:-4}" ./run
