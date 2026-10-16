#!/usr/bin/env bash
# Build ablsteady on NLR Kestrel.
#
#   tools/kestrel/build_kestrel.sh cpu   # on a CPU login node
#   tools/kestrel/build_kestrel.sh gpu   # on a GPU login node (kestrel-gpu.hpc.nrel.gov)
#
# Modules: tools/kestrel/kestrel_env.sh (ERF's Kestrel profiles at the pinned
# commit). CUDA compiles need memory; on a shared node, build inside an
# allocation with enough of it, e.g.
#   salloc -A <account> --gpus=1 --mem=240G -t 01:00:00
#
# Not yet run on Kestrel.
set -euo pipefail

target=${1:-}
case "$target" in
    cpu|gpu) ;;
    *) echo "usage: $0 cpu|gpu" >&2; exit 2 ;;
esac

here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../.." && pwd)
build="$src/build-kestrel-$target"

# shellcheck source=kestrel_env.sh
source "$here/kestrel_env.sh" "$target"

opts=(-DCMAKE_BUILD_TYPE=Release
      -DABLSTEADY_MPI=ON
      -DCMAKE_C_COMPILER=cc
      -DCMAKE_CXX_COMPILER=CC
      -DMPIEXEC_EXECUTABLE="$(command -v srun)"
      -DMPIEXEC_NUMPROC_FLAG=-n)
if [ "$target" = gpu ]; then
    # H100: compute capability 9.0.
    opts+=(-DABLSTEADY_GPU_BACKEND=CUDA
           -DAMReX_CUDA_ARCH=9.0
           -DCMAKE_CUDA_HOST_COMPILER=CC)
fi

git -C "$src" submodule update --init --recursive
cmake -S "$src" -B "$build" "${opts[@]}"
cmake --build "$build" -j "${ABLSTEADY_BUILD_JOBS:-16}"
echo "Built $build/ablsteady and $build/ablsteady_refavg"
