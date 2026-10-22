#!/bin/bash
# Runs the Wagenbrenner et al. (2019) comparison one case at a time and plots the
# transects. Usage (from the repository root):
#   Tests/Real_Terrain/WindNinja2019/run_all.sh BUILD_DIR WNTOOLS RUNDIR [NP]
# BUILD_DIR holds gust; WNTOOLS is a checkout of github.com/nwagenbrenner/windninja-tools;
# RUNDIR receives the terrain files, one directory per run, the samples and the figures.
# MODES (default "mc rans"; also mc_stable, alpha_v = 0.5, WindNinja's slightly stable class): mc is the mass-consistent model alone (steady.stage unset, the inflow
# profile projected with alpha_h = alpha_v = 1, WindNinja's neutral COM setting, on the same mesh);
# rans is the full RANS run. LIMITERS (default "van_leer") repeats each rans run per convection.limiter;
# MESH=100k replaces each deck's mesh by one of about 100 K cells, the paper's CFD budget (run on 4 ranks
# as the paper's 4 processors): runs <case>_100k_<mode>.
# CASES (default all but askervein_stable and bsb_wn) picks cases: askervein askervein_stable bolund_270
# bolund_239 bolund_90 bsb bsb_wn (askervein_stable: rans only; its mc runs are askervein's).
set -u
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
BIN=$(cd "$1" && pwd)/gust; WN=$2; RUN=$3; NP=${4:-8}
D=$ROOT/Tests/Real_Terrain/WindNinja2019
T=$ROOT/Tests/Real_Terrain/WindNinja2019/windninja2019.py
mkdir -p "$RUN" && RUN=$(cd "$RUN" && pwd)
[ -f "$RUN/masts.csv" ] || python3 "$T" prepare "$WN" "$RUN"
for c in ${CASES:-askervein bolund_270 bolund_239 bolund_90 bsb_wn_consistent}; do
  case $c in
    askervein*) terr="terrain.file=$ROOT/extern_erf/askervein.txt"; samp="--erf_terrain $ROOT/extern_erf/askervein.txt"; site=askervein; k=210; z0=0.03 ;;
    bolund_*)  terr="terrain.file=$RUN/bolund.xyz"; samp="--terrain $RUN/bolund.xyz"; site=bolund; k=${c#bolund_}; z0=0.01 ;;
    bsb_wn*)   # the paper-like mesh; bsb_wn.xyz: windninja2019.py smooth big_butte.asc --sigma 120 --cap 2050 --soft 150
               [ -f "$RUN/bsb_wn.xyz" ] || python3 "$T" smooth "$WN/test_cases/big_butte.asc" "$RUN/bsb_wn.xyz" \
                 --sigma 120 --cap 2050 --soft 150
               terr="terrain.file=$RUN/bsb_wn.xyz"; samp="--terrain $RUN/bsb_wn.xyz"; site=bsb; k=222; z0=0.01 ;;
    bsb)       terr="terrain.file=$RUN/bsb.xyz"; samp="--terrain $RUN/bsb.xyz"; site=bsb; k=222; z0=0.01 ;;
  esac
  # The 100 K mesh: n_cell, max_grid_size, dz0 and stretch_ratio (top unchanged).
  mesh100=""
  if [ "${MESH:-}" = 100k ]; then
    case $c in
      askervein*) mesh100="56 56 32|28|2.|1.1357" ;;
      bolund_*)   mesh100="80 40 32|20|0.5|1.1086" ;;
      bsb)        mesh100="56 60 32|28|2.|1.1941" ;;
      *)          echo "$c: no 100 K mesh (MESH=100k); skipped"; continue ;;
    esac
  fi
  tag=${MESH:+_$MESH}
  csvs=()
  runs=()
  for m in ${MODES:-mc rans}; do
    if [ "$m" = mc ] || [ "$m" = mc_stable ]; then runs+=("$m"); else for lim in ${LIMITERS:-van_leer}; do runs+=("rans_$lim"); done; fi
  done
  for r in "${runs[@]}"; do
    dir=$RUN/${c}${tag}_$r
    mkdir -p "$dir"
    if [ "$r" = mc ] || [ "$r" = mc_stable ]; then
      # The RANS deck without the solver's keys (steady, simplec, turb, forcing.gravity,
      # io.turb_plot_file: refused without a stage), alpha_v = alpha_h = 1.
      grep -vE "^(steady|simplec|turb|convection)\.|^init\.alpha_v|^forcing\.gravity|^io\.turb_plot_file" "$D/inputs_$c" > "$dir/inputs"
      # mc_stable: WindNinja's "slightly stable" class, alpha = alpha_H/alpha_V = 0.5 (src/ninja/stability.cpp),
      # the run the paper used for its Askervein COM result.
      if [ "$r" = mc_stable ]; then echo "init.alpha_v = 0.5" >> "$dir/inputs"; else echo "init.alpha_v = 1.0" >> "$dir/inputs"; fi
      extra=""
    else
      cp "$D/inputs_$c" "$dir/inputs"
      extra="convection.limiter=${r#rans_}"
      # upwind: the blend towards upwind has no effect and is refused, so its keys are dropped.
      if [ "$r" = rans_upwind ]; then
        grep -vE "^convection\.hybrid_" "$dir/inputs" > "$dir/inputs.tmp" && mv "$dir/inputs.tmp" "$dir/inputs"
      fi
    fi
    if [ -n "$mesh100" ]; then
      IFS='|' read -r nc mg dz sr <<< "$mesh100"
      sed -e "s/^mesh.n_cell .*/mesh.n_cell           = $nc/" -e "s/^mesh.max_grid_size .*/mesh.max_grid_size    = $mg/" \
        -e "s/^mesh.dz0 .*/mesh.dz0              = $dz/" -e "s/^mesh.stretch_ratio .*/mesh.stretch_ratio    = $sr/" \
        "$dir/inputs" > "$dir/inputs.tmp" && mv "$dir/inputs.tmp" "$dir/inputs"
    fi
    if [ ! -f "$dir/done" ]; then
      ( cd "$dir" && /usr/bin/time -l mpirun -np "$NP" "$BIN" inputs $terr $extra \
          < /dev/null > log.txt 2> time.txt ) && touch "$dir/done"
      grep -E "steady after|not steady|diverged|wrote plt" "$dir/log.txt" | head -1
    fi
    if [ -f "$dir/done" ]; then
      python3 "$T" sample "$dir/plt" --site $site --case $k --masts "$RUN/masts.csv" $samp --z0 $z0 \
        --out "$RUN/${c}${tag}_$r.csv" && csvs+=("$RUN/${c}${tag}_$r.csv")
    fi
  done
  # The plot takes every sample of the case in RUN (the mc and rans runs of earlier calls too).
  pc=${c%_stable}  # askervein_stable's samples go on the askervein plot
  all=("$RUN/${pc}"_mc.csv "$RUN/${pc}"_mc_stable.csv "$RUN/${pc}"_rans_*.csv "$RUN/${pc}"_stable_rans_*.csv
       "$RUN/${pc}"_100k_*.csv "$RUN/${pc}"_stable_100k_*.csv)
  have=()
  for f in "${all[@]}"; do [ -f "$f" ] && have+=("$f"); done
  [ ${#have[@]} -gt 0 ] && python3 "$T" plot "${have[@]}" --masts "$RUN/masts.csv" --site $site --case $k \
    --paper "$D/wagenbrenner2019_digitized.csv" --out "$RUN/$pc.png"
done
