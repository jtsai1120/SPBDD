#!/bin/bash
# Precise timing: repeat the fast tools so sub-second runs are measurable.
#   $1 prefix   $2 d   $3 timeout   $4 repeats for m4riCC
set -u
P=$1; D=$2; T=${3:-1200}; R=${4:-20}
cd /tmp/css
N=$(head -1 "$P.txt" | cut -d, -f1)

spbdd=$( { /usr/bin/time -f "@@ %e %M" timeout "$T" /tmp/dist_one "$P.txt" "$N" ; } 2>&1 )
sd=$(echo "$spbdd" | grep -oE "spbdd d=[0-9-]+" | cut -d= -f2)
st=$(echo "$spbdd" | grep "^@@" | awk '{print $2}')
sm=$(echo "$spbdd" | grep "^@@" | awk '{printf "%.0f", $3/1024}')
[ -z "$sd" ] && { sd="timeout"; st=">$T"; sm="-"; }

# m4riCC: exact connected cluster, run R times so the timer resolves it
mt=$( { /usr/bin/time -f "@@ %e %M" bash -c \
        "for i in \$(seq $R); do /tmp/dist-m4ri/src/dist_m4ri method=2 wmax=$D \
           finH=${P}_HX.mtx finG=${P}_HZ.mtx >/dev/null 2>&1; done" ; } 2>&1 | grep "^@@")
one=$(echo "$mt" | awk -v r="$R" '{printf "%.4f", $2/r}')
mm=$(echo "$mt" | awk '{printf "%.0f", $3/1024}')
md=$(/tmp/dist-m4ri/src/dist_m4ri method=2 wmax="$D" "finH=${P}_HX.mtx" "finG=${P}_HZ.mtx" 2>/dev/null | tail -1 | awk '{print $1}')

printf "%-12s n=%-4s d=%-3s | SPBDD %10s s  d=%-8s %4s MB | m4riCC %9s s  d=%-4s %3s MB\n" \
       "$P" "$N" "$D" "$st" "$sd" "$sm" "$one" "${md:-?}" "$mm"
