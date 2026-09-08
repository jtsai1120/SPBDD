# Distance finding: SPBDD against the tools from arXiv:2603.22532

All timings were taken on one machine (i5-12600KF, 16 threads, 16 GB, WSL2
Ubuntu 22.04), single-threaded except `m4riCC`, which is multithreaded C. The
paper does not disclose the hardware or the timeouts it used anywhere in its
text, so its published tables are not a sound basis for comparison; everything
below was re-measured here.

Every tool was given byte-identical input: the generators in `gen_*.py` write
both the codetables.de symplectic form (for SPBDD and the Python methods) and
the MatrixMarket pair (for `dist_m4ri`).

`build_tools.sh` builds the reference tools without root.

## What the tools are

| tool | method | exact? | takes non-CSS? |
|---|---|---|---|
| `BZDistMW` | Brouwer-Zimmermann | yes | yes |
| `connectedClusterMW` | connected cluster, pure Python | yes | yes |
| `pySATDist` | MaxSAT | yes | yes |
| `QDistEvol` | evolutionary search | no, upper bound | yes |
| `m4riCC` | connected cluster, multithreaded C | yes | **no, CSS only** |
| `m4riRW` | random window | no, upper bound | **no, CSS only** |

`dist-m4ri` is what the paper recommends for quantum CSS codes, and it cannot
read a non-CSS code at all: its `css` flag is documented as "reserved for
future use (1)" and it computes `min(d_X, d_Z)` from the two classical halves.

## Non-CSS: codetables.de best-known codes

Dense, unstructured -- the worst case for a decision diagram.

| code | d | SPBDD | BZDistMW | connCluster | pySATDist | QDistEvol |
|---|---|---|---|---|---|---|
| `[[19,5,5]]` | 5 | 0.005 | 0.021 | 0.207 | 0.119 | 0.873 |
| `[[21,3,6]]` | 6 | 0.020 | 0.062 | 3.976 | 0.458 | 0.996 |
| `[[27,6,6]]` | 6 | **0.060** | 0.454 | -- | 2.565 | 0.945 |
| `[[29,6,7]]` | 7 | **0.681** | 2.968 | -- | 17.260 | 1.053 |
| `[[30,9,6]]` | 6 | **0.359** | 1.387 | -- | 4.830 | 1.128 |
| `[[32,5,8]]` | 8 | **4.896** | 17.908 | -- | timeout | 1.113 |
| `[[33,7,8]]` | 8 | **12.537** | 35.573 | -- | timeout | 1.365 |
| `[[34,17,5]]` | 5 | **0.152** | 50.270 | -- | -- | 1.202 |
| `[[40,5,9]]` | 9 | **316.7** | 453.0 | -- | -- | 1.749 |
| `[[41,2,11]]` | 11 | timeout >3600 | 451.6 | -- | -- | 1.757 |
| `[[49,17,7]]` | 7 | **60.0** | 359.3 | -- | -- | 1.641 |

SPBDD leads over a window roughly `n = 27..49` while `d <= 9`, and falls off a
cliff at `d >= 11`. `QDistEvol` is near-constant but only ever returns an upper
bound; it happened to be right on every code here.

## CSS: toric, bivariate bicycle, quantum Tanner

Sparse and local -- and the case `m4riCC` is built for.

| code | n | d | SPBDD | m4riCC | BZDistMW (CSS path) |
|---|---|---|---|---|---|
| Tanner `[[36,8,4]]` | 36 | 4 | 0.03 | **0.0038** | -- |
| toric `[[32,2,4]]` | 32 | 4 | 0.02 | **0.0035** | 0.002 |
| toric `[[50,2,5]]` | 50 | 5 | 0.14 | **0.0040** | 0.004 |
| toric `[[72,2,6]]` | 72 | 6 | 1.93 | **0.0045** | 0.029 |
| BB `[[72,12,6]]` | 72 | 6 | 3.13 | **0.0045** | -- |
| toric `[[98,2,7]]` | 98 | 7 | 25.0 | **0.0045** | 20.4 |
| toric `[[128,2,8]]` | 128 | 8 | 319.5 | **<0.001** | timeout >600 |
| BB `[[108,8,10]]` | 108 | 10 | timeout >1800 | **<0.001** | -- |

`m4riCC` is flat at about 4 ms from n=32 to n=128 and d=4 to d=10, on 7-11 MB.
SPBDD grows by roughly 12x per unit of distance, on 32-49 MB.

Two fairness checks were run, because the first version of each comparison was
wrong:

- BZDistMW was first given the toric codes through the non-CSS path, where it
  took 77 s on `[[72,2,6]]`. Through the CSS path -- which is what one would
  actually use on a CSS code -- it takes 0.029 s. The 2600x difference was an
  artefact of the harness, not a property of the tool.
- `m4riCC` was first given `wmax` equal to the answer. Raising it to 20 on
  `[[72,2,6]]` changes nothing (0.0040 s either way), because CC stops as soon
  as it has certified the minimum. The hint is worth nothing and the timings
  stand.

## Why the gap

A decision diagram for `{e : wt(e) <= d} & N(S)` has to distinguish, at every
level, every reachable pair of (weight so far, partial syndrome). The number of
those that survive the weight budget grows like `C(n, d)`, which is why the
cost is exponential in the distance rather than in n -- `[[49,17,7]]` finishes
in a minute while the smaller `[[41,2,11]]` does not finish in an hour.

Connected cluster never forms that set. A minimum-weight logical operator of a
sparse code is connected in the Tanner graph -- were it not, each component
would be a lighter logical operator -- so CC only enumerates *connected*
supports, which on a sparse graph is a vanishing fraction of the `C(n, d)`
supports the diagram represents. It exploits the same locality a diagram hopes
to, and far more directly.

## Where this leaves SPBDD

Distance finding is not the case for a BDD representation. The one position
that holds up is **non-CSS codes at small d**: the tools that dominate the CSS
benchmark cannot read a non-CSS code at all, and SPBDD treats CSS and non-CSS
identically because it works on the symplectic form throughout. That window is
narrow, but the measurements above are inside it.

## Reproducing

```bash
bash experiments/build_tools.sh          # pip + codeDistance + m4ri + dist-m4ri, no root
make && g++ -std=c++17 -O2 -Iinclude experiments/distance_one.cpp \
    build/libspbdd.a -lm -o /tmp/dist_one

python3 experiments/gen_toric.py 2 8 > /tmp/toric.txt
python3 experiments/split_css.py /tmp/toric.txt /tmp/css
python3 experiments/gen_css.py bb 6 6 /tmp/css/bb72 6

bash experiments/css_bench.sh toric98 7 1200 20
python3 experiments/run_tools.py QECC_sample.txt 5 26 60 BZDistMW QDistEvol
python3 experiments/run_css.py  toric.txt 8 128 600 BZDistMW
```
