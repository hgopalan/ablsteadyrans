#!/usr/bin/env bash
# Regenerate Tests/Gold/erf_hill2d_{btf,stf,sullivan} from ERF itself.
#
#   tools/make_erf_mesh_gold.sh /path/to/erf_exec /path/to/ERF
#
# ERF must be checked out at the commit in external/VERSIONS and built
# (MPI and FFT on). Each run is ERF's Neutral_Hill_2D deck shrunk to
# 32 x 1 x 16 cells and stopped at step 0, so only the mesh is computed.
set -euo pipefail
exe=$1; erf=$2
src=$(cd "$(dirname "$0")/.." && pwd)
want=$(awk '$1 == "erf" {print $2}' "$src/external/VERSIONS")
have=$(git -C "$erf" rev-parse HEAD)
[ "$have" = "$want" ] || { echo "ERF is at $have, external/VERSIONS pins $want" >&2; exit 1; }
work=$(mktemp -d)
for n in btf:0 stf:1 sullivan:2; do
    name=${n%%:*}; sm=${n##*:}
    mkdir -p "$work/$name"
    cp "$erf/Exec/CanonicalTests/Canonical_RANS/Neutral_Hill_2D"/* "$work/$name/"
    (cd "$work/$name" && "$exe" inputs_hill max_step=0 "amr.n_cell=32 1 16" erf.terrain_smoothing=$sm \
        erf.plot_int_1=1 "erf.plot_vars_1=z_phys detJ" erf.check_int=-1 erf.v=0 > run.log)
    rm -rf "$src/Tests/Gold/erf_hill2d_$name"
    mkdir -p "$src/Tests/Gold/erf_hill2d_$name"
    cp -R "$work/$name/plt00000/Header" "$work/$name/plt00000/Level_0" "$src/Tests/Gold/erf_hill2d_$name/"
done
echo "Wrote Tests/Gold/erf_hill2d_{btf,stf,sullivan} (work files in $work)"
