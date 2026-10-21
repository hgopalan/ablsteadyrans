#!/usr/bin/env bash
# Dry run of tools/kestrel/erf_references.sh with stub module/sbatch/srun.
#
#   kestrel_dryrun.sh SRC REFAVG PLOTFILE WORK
#
# submit: the Askervein job script carries the run length, the plotfile
#   period, the 2-D plotfile switch-off and the rank count.
# archive: of a fake finished run holding 3-D plotfiles inside and outside
#   the final window and 2-D plotfiles under a similar prefix, only the 3-D
#   ones in the window are averaged, and the archive verifies.
set -euo pipefail
src=$1; refavg=$2; plt=$3; work=$4

rm -rf "$work"
mkdir -p "$work/bin" "$work/ERF/Exec/CanonicalTests/Real_Terrain/Askervein"
for stub in module sbatch srun; do
    printf '#!/bin/sh\necho "%s stub: $*"\n' "$stub" > "$work/bin/$stub"
    chmod +x "$work/bin/$stub"
done
case_dir="$work/ERF/Exec/CanonicalTests/Real_Terrain/Askervein"
printf 'erf.plot_file_1     = "plt_"   # 3-D\nerf.plot2d_file_1 = "plt2d_"\n' > "$case_dir/inputs_anel"
echo "terrain" > "$case_dir/askervein.txt"

export PATH="$work/bin:$PATH" GUST_REF_WORK="$work" GUST_ACCOUNT=test GUST_REFAVG="$refavg"
script="$src/tools/kestrel/erf_references.sh"

"$script" submit cpu Askervein
run="$work/runs/Askervein_cpu"
grep -q 'srun ".*/erf_exec" "inputs_anel" stop_time=4800 erf.plot_per_1=100 erf.plot_int_1=-1 erf.plot2d_int_1=-1' "$run/job.sh" \
    || { echo "job.sh lacks the expected overrides:"; cat "$run/job.sh"; exit 1; }
grep -q '^#SBATCH --ntasks=104$' "$run/job.sh" || { echo "job.sh lacks --ntasks=104"; exit 1; }
[ "$(cat "$run/inputs_list" | tr '\n' ' ')" = "askervein.txt inputs_anel overrides " ] \
    || { echo "unexpected inputs_list: $(cat "$run/inputs_list")"; exit 1; }

# The finished run: plt_ snapshots at 1000 s (outside the window) and 2400 to
# 4800 s, 2-D plotfiles plt2d_ (other variables), and the run's bookkeeping.
copy_at () {
    cp -R "$plt" "$run/case/$1"
    awk -v t="$2" 'NR == 2 { n = $1 } NR == 4 + n { print t; next } { print }' "$plt/Header" > "$run/case/$1/Header"
}
copy_at plt_00100 1000
i=0
for t in 2400 3000 3600 4200 4800; do i=$((i + 1)); copy_at "plt_0$((1000 + i))" "$t"; done
copy_at plt2d_00010 4500
sed -i.bak '3s/.*/u_star/' "$run/case/plt2d_00010/Header" && rm -f "$run/case/plt2d_00010/Header.bak"
echo "surface history" > "$run/case/surf_hist.dat"   # a run output, not an input
echo 1234 > "$run/wall_seconds"

"$script" archive cpu Askervein | tee "$work/archive.log"
grep -q "verifies: case Askervein, ERF 92410d578, 5 snapshots, kestrel-cpu, 1 node(s), 104 rank(s)" "$work/archive.log" \
    || { echo "archive did not average exactly the 5 plotfiles in the window"; exit 1; }
for f in inputs/inputs_anel inputs/askervein.txt inputs/overrides; do
    [ -f "$work/archive/Askervein_cpu/$f" ] || { echo "archive lacks $f"; exit 1; }
done
[ ! -e "$work/archive/Askervein_cpu/inputs/surf_hist.dat" ] || { echo "a run output was archived as an input"; exit 1; }
# A deck that keeps the default prefix "plt": the 2-D plotfiles plt2d_* match
# the glob plt*, and only the name pattern keeps them out.
flat="$work/ERF/Exec/CanonicalTests/Canonical_RANS/Neutral_ABL_Flat"
mkdir -p "$flat"
echo "erf.plot_vars_1 = x_velocity" > "$flat/inputs_neutral"
"$script" submit cpu Neutral_ABL_Flat
run="$work/runs/Neutral_ABL_Flat_cpu"
i=0
for t in 169200 171000 172800; do i=$((i + 1)); copy_at "plt0$((1000 + i))" "$t"; done
copy_at plt2d_00020 172800
sed -i.bak '3s/.*/u_star/' "$run/case/plt2d_00020/Header" && rm -f "$run/case/plt2d_00020/Header.bak"
echo 600 > "$run/wall_seconds"
"$script" archive cpu Neutral_ABL_Flat | tee "$work/archive_flat.log"
grep -q "verifies: case Neutral_ABL_Flat, ERF 92410d578, 3 snapshots" "$work/archive_flat.log" \
    || { echo "the default-prefix archive did not average exactly the 3 plotfiles"; exit 1; }

echo "kestrel dry run passed"
