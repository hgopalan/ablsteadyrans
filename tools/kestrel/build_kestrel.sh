#!/usr/bin/env bash
# Build ablsteady on NLR Kestrel.
#
#   tools/kestrel/build_kestrel.sh cpu   # on a CPU login node
#   tools/kestrel/build_kestrel.sh gpu   # on a GPU login node (kestrel-gpu.hpc.nrel.gov)
#
# Modules follow ERF's Kestrel guide (Docs/sphinx_doc/kestrel_build_run.rst at
# ERF 92410d5). Kestrel's module names change with system updates: check them
# with `module avail` before building and update this script if they moved.
#
# CUDA compiles need memory. On a shared node, build inside an allocation with
# enough of it, e.g.
#   salloc -A <account> -p debug --gres=gpu:1 --mem=240G -t 01:00:00
#
# Not yet run on Kestrel: the first run is part of sub-phase 1.1 (PR-03).
set -euo pipefail

target=${1:-}
case "$target" in
    cpu|gpu) ;;
    *) echo "usage: $0 cpu|gpu" >&2; exit 2 ;;
esac

src=$(cd "$(dirname "$0")/../.." && pwd)
build="$src/build-kestrel-$target"

module purge
module load PrgEnv-gnu/8.5.0
module load cmake

opts=(-DCMAKE_BUILD_TYPE=Release
      -DABLSTEADY_MPI=ON
      -DCMAKE_C_COMPILER=cc
      -DCMAKE_CXX_COMPILER=CC
      -DMPIEXEC_EXECUTABLE="$(command -v srun)"
      -DMPIEXEC_NUMPROC_FLAG=-n)

if [ "$target" = gpu ]; then
    module load cuda/12.3
    module load craype-x86-milan
    # H100: compute capability 9.0.
    opts+=(-DABLSTEADY_GPU_BACKEND=CUDA
           -DAMReX_CUDA_ARCH=9.0
           -DCMAKE_CUDA_HOST_COMPILER=CC)
fi

git -C "$src" submodule update --init --recursive
cmake -S "$src" -B "$build" "${opts[@]}"
cmake --build "$build" -j "${ABLSTEADY_BUILD_JOBS:-16}"
echo "Built $build/ablsteady"
