#!/bin/bash
# One-change sensitivity runs of a WindNinja2019 case (GUST - CFD), then the samples, a plot and the statistics.
# Usage (from the repository root, after run_all.sh has prepared RUNDIR):
#   Tests/Real_Terrain/WindNinja2019/sweep.sh BUILD_DIR RUNDIR CASE VARIANTS [NP]
# VARIANTS is a file of lines "name key=value [key=value ...]"; a key the deck sets is replaced in the
# run's copy of the deck, any other key is added, and key=- removes it. Each run goes to RUNDIR/<CASE>_rans__<name>; its
# sample is RUNDIR/<CASE>_rans__<name>.csv; the plot RUNDIR/<CASE>_sweep.png shows the baseline
# (<CASE>_rans_van_leer.csv) and every variant; the statistics go to RUNDIR/<CASE>_sweep_stats.md.
set -u
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
BIN=$(cd "$1" && pwd)/gust; RUN=$(cd "$2" && pwd); c=$3; VARS=$4; NP=${5:-8}
D=$ROOT/Tests/Real_Terrain/WindNinja2019
T=$ROOT/Tests/Real_Terrain/WindNinja2019/windninja2019.py
case $c in
  askervein*) terr="terrain.file=$ROOT/extern_erf/askervein.txt"; samp="--erf_terrain $ROOT/extern_erf/askervein.txt"; site=askervein; k=210; z0=0.03 ;;
  bolund_*)   terr="terrain.file=$RUN/bolund.xyz"; samp="--terrain $RUN/bolund.xyz"; site=bolund; k=${c#bolund_}; z0=0.01 ;;
  bsb_wn*)    terr="terrain.file=$RUN/bsb_wn.xyz"; samp="--terrain $RUN/bsb_wn.xyz"; site=bsb; k=222; z0=0.01 ;;
  bsb|bsb_srtm) terr="terrain.file=$RUN/bsb.xyz"; samp="--terrain $RUN/bsb.xyz"; site=bsb; k=222; z0=0.01 ;;
esac
csvs=("$RUN/${c}_rans_van_leer.csv")
while read -r name kv; do
  [ -z "$name" ] || [ "${name:0:1}" = "#" ] && continue
  dir=$RUN/${c}_rans__$name
  mkdir -p "$dir"
  cp "$D/inputs_$c" "$dir/inputs"
  zs=$z0
  for p in $kv; do
    key=${p%%=*}; val=${p#*=}
    if [ "$val" = - ]; then
      grep -v "^$key " "$dir/inputs" > "$dir/inputs.tmp" && mv "$dir/inputs.tmp" "$dir/inputs"
    elif grep -q "^$key " "$dir/inputs"; then
      sed "s|^$key .*|$key = $val|" "$dir/inputs" > "$dir/inputs.tmp" && mv "$dir/inputs.tmp" "$dir/inputs"
    else
      echo "$key = $val" >> "$dir/inputs"
    fi
    [ "$key" = bc.surface_z0 ] && zs=$val
  done
  if [ ! -f "$dir/done" ]; then
    ( cd "$dir" && /usr/bin/time -l mpirun -np "$NP" "$BIN" inputs $terr < /dev/null > log.txt 2> time.txt ) && touch "$dir/done"
  fi
  echo "$name: $(grep -E 'steady after|not steady|diverged' "$dir/log.txt" | head -1) $(grep ' real' "$dir/time.txt" | awk '{print $1 " s"}')"
  [ -f "$dir/done" ] && python3 "$T" sample "$dir/plt" --site $site --case $k --masts "$RUN/masts.csv" $samp --z0 "$zs" \
      --out "$RUN/${c}_rans__$name.csv" > /dev/null && csvs+=("$RUN/${c}_rans__$name.csv")
done < "$VARS"
python3 "$T" plot "${csvs[@]}" --masts "$RUN/masts.csv" --site $site --case $k \
  --paper "$D/wagenbrenner2019_digitized.csv" --out "$RUN/${c}_sweep.png"
python3 "$T" stats "${csvs[@]}" --masts "$RUN/masts.csv" --site $site > "$RUN/${c}_sweep_stats.md"
