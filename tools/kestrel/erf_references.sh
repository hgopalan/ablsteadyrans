#!/usr/bin/env bash
# ERF reference solutions on NLR Kestrel.
#
#   erf_references.sh build   cpu|gpu          build ERF at the pinned commit
#   erf_references.sh smoke   cpu|gpu          ERF's `ctest -L rans`
#                                              (run inside an allocation)
#   erf_references.sh submit  cpu|gpu CASE...  submit reference runs ("all" = every case)
#   erf_references.sh archive cpu|gpu CASE     time-average the finished run into an
#                                              archive (run on a CPU login node)
#
# Environment:
#   GUST_ACCOUNT   Kestrel allocation (required for submit)
#   GUST_REF_WORK  work directory (default /scratch/$USER/gust-refs)
#   GUST_REFAVG    gust_refavg to use (default build-kestrel-cpu/gust_refavg)
#
# Each run writes a plotfile every window/24 seconds of physical time
# (erf.plot_per_1); the archive is the trapezoid time average of the
# plotfiles in the final window, with the srun wall time, node and rank count,
# the deck, its input files and the command-line overrides.
#
# A run that hits its time limit is rerun from scratch (remove its run
# directory); there is no restart support.
#
# Not yet run on Kestrel. Partition and GPU request syntax follow Kestrel's
# documentation as of October 2026 and may need adjusting.
set -euo pipefail

# Called ABLSTEADY_* before the rename to GUST (October 2026); the old names still work.
GUST_ACCOUNT=${GUST_ACCOUNT:-${ABLSTEADY_ACCOUNT:-}}
[ -n "$GUST_ACCOUNT" ] || unset GUST_ACCOUNT

here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../.." && pwd)
work=${GUST_REF_WORK:-${ABLSTEADY_REF_WORK:-/scratch/$USER/gust-refs}}
erf_sha=$(awk '$1 == "erf" {print $2}' "$src/external/VERSIONS")
erf_src="$work/ERF"
refavg=${GUST_REFAVG:-${ABLSTEADY_REFAVG:-$src/build-kestrel-cpu/gust_refavg}}

usage () { sed -n '2,25p' "$0" >&2; exit 2; }

# CASE -> ERF directory (under Exec/CanonicalTests), deck, stop time (s),
# averaging window (s), MPI ranks (CPU) and GPUs, wall-clock limit.
#
# The neutral flat case and both hills start from a uniform geostrophic wind
# with Coriolis on, which sets off an inertial oscillation (period 2 pi / f =
# 17.5 h at f = 1e-4 s^-1). They run 48 h (2.75 periods) so the final window
# is close to steady, instead of the decks' 12, 6 and 4 h. The stable (GABLS1,
# surface cooling) and convective (growing mixed layer) cases are transient
# by design: they keep their decks' lengths and are not steady references.
case_info () {
    case "$1" in
        Neutral_ABL_Flat)    echo "Canonical_RANS/Neutral_ABL_Flat inputs_neutral 172800 3600 8 1 06:00:00" ;;
        Stable_ABL_Flat)     echo "Canonical_RANS/Stable_ABL_Flat inputs_stable 32400 3600 8 1 02:00:00" ;;
        Convective_ABL_Flat) echo "Canonical_RANS/Convective_ABL_Flat inputs_convective 14400 3600 8 1 02:00:00" ;;
        Neutral_Hill_2D)     echo "Canonical_RANS/Neutral_Hill_2D inputs_hill 172800 3600 32 1 16:00:00" ;;
        Neutral_Hill_3D)     echo "Canonical_RANS/Neutral_Hill_3D inputs_hill3d 172800 3600 64 1 24:00:00" ;;
        # The deck stops at 600 s; the lee flow needs about 2000 s to pass
        # through. The reference runs 4800 s and averages the last 2400 s.
        Askervein)           echo "Real_Terrain/Askervein inputs_anel 4800 2400 104 4 12:00:00" ;;
        *) echo "unknown case: $1" >&2; exit 2 ;;
    esac
}
all_cases="Neutral_ABL_Flat Stable_ABL_Flat Convective_ABL_Flat Neutral_Hill_2D Neutral_Hill_3D Askervein"

# The plotfile prefix a deck sets (erf.plot_file_1, default plt).
plot_prefix () {
    local p
    p=$(awk -F= '$1 ~ /^[ \t]*erf\.plot_file_1[ \t]*$/ { gsub(/[ \t"]/, "", $2); sub(/#.*/, "", $2); print $2 }' "$1" | tail -n 1)
    echo "${p:-plt}"
}

cmd=${1:-}; target=${2:-}
[ -n "$cmd" ] || usage
case "$target" in cpu|gpu) ;; *) usage ;; esac
shift 2

erf_build="$work/build-erf-$target"
erf_exe="$erf_build/Exec/erf_exec"

case "$cmd" in
build)
    # shellcheck source=kestrel_env.sh
    source "$here/kestrel_env.sh" "$target"
    mkdir -p "$work"
    if [ ! -d "$erf_src/.git" ]; then
        git clone https://github.com/erf-model/ERF.git "$erf_src"
    fi
    git -C "$erf_src" fetch origin
    git -C "$erf_src" checkout --detach "$erf_sha"
    git -C "$erf_src" submodule update --init Submodules/AMReX
    opts=(-DCMAKE_BUILD_TYPE=Release
          -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC
          -DERF_ENABLE_MPI=ON -DERF_ENABLE_FFT=ON
          -DERF_ENABLE_TESTS=ON -DERF_ENABLE_UNIT_TESTS=OFF -DERF_ENABLE_FCOMPARE=ON
          -DMPIEXEC_EXECUTABLE="$(command -v srun)" -DMPIEXEC_NUMPROC_FLAG=-n)
    if [ "$target" = gpu ]; then
        opts+=(-DERF_ENABLE_CUDA=ON -DAMReX_CUDA_ARCH=9.0 -DCMAKE_CUDA_HOST_COMPILER=CC)
    fi
    cmake -S "$erf_src" -B "$erf_build" "${opts[@]}"
    cmake --build "$erf_build" -j "${GUST_BUILD_JOBS:-${ABLSTEADY_BUILD_JOBS:-16}}"
    echo "ERF $erf_sha built in $erf_build"
    ;;

smoke)
    # Guards the reference generation. Needs compute nodes.
    if [ -z "${SLURM_JOB_ID:-}" ]; then
        echo "run inside an allocation, e.g. salloc -A \$GUST_ACCOUNT -N 1 -t 01:00:00" >&2
        exit 2
    fi
    # shellcheck source=kestrel_env.sh
    source "$here/kestrel_env.sh" "$target"
    ctest --test-dir "$erf_build" -L rans --output-on-failure
    ;;

submit)
    : "${GUST_ACCOUNT:?set GUST_ACCOUNT to the Kestrel allocation}"
    [ "$#" -gt 0 ] || usage
    [ "$1" = all ] && set -- $all_cases
    for c in "$@"; do
        read -r dir deck stop window ranks gpus walltime <<< "$(case_info "$c")"
        run="$work/runs/${c}_$target"
        if [ -e "$run" ]; then echo "$run exists; remove it to rerun $c" >&2; exit 1; fi
        mkdir -p "$run/case"
        cp "$erf_src/Exec/CanonicalTests/$dir"/* "$run/case"/
        per=$(awk -v w="$window" 'BEGIN { printf "%.10g", w / 24 }')
        # 2-D plotfiles off: Askervein writes them every 10 steps under a
        # prefix the 3-D plotfiles share.
        overrides="stop_time=$stop erf.plot_per_1=$per erf.plot_int_1=-1 erf.plot2d_int_1=-1"
        if [ "$target" = gpu ]; then
            n=$gpus
            res="#SBATCH --ntasks=$gpus
#SBATCH --gpus=$gpus"
        else
            n=$ranks
            res="#SBATCH --ntasks=$ranks"
        fi
        {
            echo "$erf_sha" > "$run/erf_sha"
            echo 1 > "$run/nodes"
            echo "$n" > "$run/ranks"
            echo "$stop $window" > "$run/window"
            echo "$overrides" > "$run/case/overrides"
            # The inputs as submitted; the run adds its own output files to case/.
            (cd "$run/case" && ls -1) > "$run/inputs_list"
        }
        cat > "$run/job.sh" <<JOB
#!/usr/bin/env bash
#SBATCH --account=$GUST_ACCOUNT
#SBATCH --job-name=ref-$c-$target
#SBATCH --time=$walltime
#SBATCH --nodes=1
$res
#SBATCH --output=$run/run.log
set -euo pipefail
source "$here/kestrel_env.sh" $target
cd "$run/case"
start=\$(date +%s)
srun "$erf_exe" "$deck" $overrides
echo \$(( \$(date +%s) - start )) > "$run/wall_seconds"
JOB
        (cd "$run" && sbatch job.sh)
    done
    ;;

archive)
    [ "$#" -eq 1 ] || usage
    c=$1
    read -r _ deck _ _ _ _ _ <<< "$(case_info "$c")"
    run="$work/runs/${c}_$target"
    for f in wall_seconds window nodes ranks erf_sha inputs_list; do
        [ -f "$run/$f" ] || { echo "$run/$f missing: the run did not finish" >&2; exit 1; }
    done
    read -r stop window < "$run/window"
    prefix=$(plot_prefix "$run/case/$deck")
    # 3-D plotfiles of the deck's prefix whose time (Header line 4 + ncomp)
    # lies in the final window.
    plts=$(cd "$run/case" && for p in "$prefix"*; do
        [[ "$p" =~ ^${prefix}[0-9]+$ ]] || continue
        [ -f "$p/Header" ] || continue
        t=$(awk 'NR == 2 { n = $1 } NR == 4 + n { print; exit }' "$p/Header")
        awk -v t="$t" -v s="$stop" -v w="$window" 'BEGIN { exit !(t >= s - w - 1e-6 * w) }' && echo "$p"
    done | sort)
    [ -n "$plts" ] || { echo "no ${prefix}* plotfiles in the final window of $run" >&2; exit 1; }
    inputs=$(tr '\n' ' ' < "$run/inputs_list")
    out="$work/archive/${c}_$target"
    mkdir -p "$work/archive"
    # The archive tool is a CPU build and runs under the CPU environment,
    # also for GPU runs (the GPU environment's MPICH_GPU_SUPPORT_ENABLED
    # would stop a binary built without GPU support).
    (
        # shellcheck source=kestrel_env.sh
        source "$here/kestrel_env.sh" cpu
        cd "$run/case"
        "$refavg" "refavg.plotfiles=$(echo $plts)" "refavg.output=$out" "refavg.case=$c" \
            "refavg.erf_sha=$(cat "$run/erf_sha")" "refavg.machine=kestrel-$target" \
            "refavg.nodes=$(cat "$run/nodes")" "refavg.ranks=$(cat "$run/ranks")" \
            "refavg.wall_seconds=$(cat "$run/wall_seconds")" "refavg.inputs_files=$inputs"
        "$refavg" "refavg.verify=$out"
    )
    ;;

*)
    usage
    ;;
esac
