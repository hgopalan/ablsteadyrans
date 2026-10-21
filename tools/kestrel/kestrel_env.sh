#!/usr/bin/env bash
# Kestrel modules for GUST and for the ERF reference builds:
#
#   source tools/kestrel/kestrel_env.sh cpu   # CPU nodes (Sapphire Rapids)
#   source tools/kestrel/kestrel_env.sh gpu   # GPU nodes (H100 + Genoa)
#
# Taken from ERF's Build/machines/kestrel_erf_cpu.profile and
# kestrel_erf_gpu.profile at the pinned ERF commit (external/VERSIONS), which
# follow Kestrel's RHEL9 upgrade (erf-model/ERF#4292). The HDF5 and NetCDF
# modules of those profiles are left out: neither code is built with them
# here. Module names change with system updates; check with `module avail`.
#
# Not yet run on Kestrel.

case "${1:-}" in
    cpu)
        module reset
        module load cpe-stack
        module load PrgEnv-gnu
        module load craype-x86-spr
        module load cmake
        module load cray-fftw/3.3.10.10
        ;;
    gpu)
        module reset
        module load cpe-stack
        module load PrgEnv-gnu
        module load craype-x86-genoa
        module load cmake
        module load cudatoolkit
        module load craype-accel-nvidia90
        module load cray-fftw/3.3.10.10
        export MPICH_GPU_SUPPORT_ENABLED=1
        ;;
    *)
        echo "usage: source kestrel_env.sh cpu|gpu" >&2
        return 2 2>/dev/null || exit 2
        ;;
esac
