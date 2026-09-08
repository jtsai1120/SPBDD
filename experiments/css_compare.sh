#!/bin/bash
# One CSS code, every tool that can take it, on this machine.
#
#   $1  prefix under /tmp/css   (expects <p>.txt, <p>_HX.mtx, <p>_HZ.mtx)
#   $2  the code's distance, so m4riCC can be given a wmax
#   $3  per-tool timeout in seconds
set -u
P=$1; D=$2; T=${3:-600}
cd /tmp/css

echo "### $P  (d=$D) ###"

# --- SPBDD ------------------------------------------------------------------
S=$( { /usr/bin/time -f "%e %M" timeout "$T" /tmp/dist_one "$P.txt" \
        "$(head -1 "$P.txt" | cut -d, -f1)" ; } 2>&1 )
echo "SPBDD              : $(echo "$S" | head -1)"
echo "                     $(echo "$S" | tail -1 | awk '{printf "%.1f s wall, %.0f MB peak RSS", $1, $2/1024}')"

# --- m4riCC : connected cluster, exact, the paper's pick for CSS -------------
# d = min(dX, dZ), so both orientations are run.
for pair in "HX HZ" "HZ HX"; do
    set -- $pair
    R=$( { /usr/bin/time -f "%e %M" timeout "$T" /tmp/dist-m4ri/src/dist_m4ri \
             method=2 wmax="$D" "finH=${P}_$1.mtx" "finG=${P}_$2.mtx" ; } 2>&1 )
    OUT=$(echo "$R" | grep -E "^[0-9]+ [0-9]+ [0-9]+$" | tail -1)
    TIME=$(echo "$R" | tail -1)
    echo "m4riCC ($1 side)  : ${OUT:-no result}   $(echo "$TIME" | awk '{printf "%.1f s wall, %.0f MB peak RSS", $1, $2/1024}')"
done
echo
