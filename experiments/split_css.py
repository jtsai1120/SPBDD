"""Split a codetables-format CSS code into HX/HZ MatrixMarket files, so
dist_m4ri and SPBDD read the same code.

  usage: split_css.py <in.txt> <out-prefix-dir>
"""
import os
import sys


def write_mtx(path, rows, nc):
    nnz = sum(sum(r) for r in rows)
    with open(path, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{len(rows)} {nc} {nnz}\n")
        for i, r in enumerate(rows):
            for j, v in enumerate(r):
                if v:
                    fh.write(f"{i + 1} {j + 1} 1\n")


src, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

params, rows = None, []
blocks = []
for line in open(src):
    line = line.strip()
    if not line:
        if params and rows:
            blocks.append((params, rows))
        params, rows = None, []
        continue
    if "," in line:
        if params and rows:
            blocks.append((params, rows))
            rows = []
        params = [int(x) for x in line.split(",")]
    elif params is not None:
        rows.append([int(c) for c in line])
if params and rows:
    blocks.append((params, rows))

for (n, k, d), rs in blocks:
    HX = [r[:n] for r in rs if any(r[:n]) and not any(r[n:])]
    HZ = [r[n:] for r in rs if any(r[n:]) and not any(r[:n])]
    if len(HX) + len(HZ) != len(rs):
        print(f"[[{n},{k},{d}]] is not CSS, skipped")
        continue
    p = os.path.join(outdir, f"toric{n}")
    with open(p + ".txt", "w") as fh:
        fh.write(f"{n},{k},{d}\n")
        for r in rs:
            fh.write("".join(map(str, r)) + "\n")
        fh.write("\n")
    write_mtx(p + "_HX.mtx", HX, n)
    write_mtx(p + "_HZ.mtx", HZ, n)
    print(f"toric{n}: [[{n},{k},{d}]]  HX {len(HX)}x{n}  HZ {len(HZ)}x{n}")
