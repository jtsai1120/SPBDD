"""Emit CSS codes in two forms from the same source, so every tool sees the
same code:

  <name>.txt          codetables.de QECC format (symplectic) -- for SPBDD and
                      the codeDistance package
  <name>_HX.mtx       MatrixMarket -- for dist_m4ri
  <name>_HZ.mtx

Sources:
  tanner  <HX.mtx> <HZ.mtx>       the quantum Tanner codes shipped with the
                                  codeDistance package
  bb      <l> <m> <A...> <B...>   bivariate bicycle, built from scratch since
                                  the package's generator needs Python 3.12

  usage: gen_css.py tanner HX.mtx HZ.mtx <out-prefix> <d>
         gen_css.py bb <l> <m> <out-prefix> <d>
"""
import sys


def read_mtx(path):
    """MatrixMarket coordinate -> dense list of rows of 0/1."""
    rows = []
    with open(path) as fh:
        header_seen = False
        for line in fh:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if not header_seen:
                nr, nc = int(parts[0]), int(parts[1])
                rows = [[0] * nc for _ in range(nr)]
                header_seen = True
                continue
            i, j = int(parts[0]) - 1, int(parts[1]) - 1
            rows[i][j] ^= 1
    return rows


def bb_matrices(l, m, a_powers, b_powers):
    """Bivariate bicycle: A = sum of x^i y^j over a_powers, likewise B.

    x is the cyclic shift on the l block, y on the m block; qubits are indexed
    (i, j) -> i * m + j, and the code has n = 2 * l * m.
    """
    nb = l * m

    def circ(powers):
        mat = [[0] * nb for _ in range(nb)]
        for r in range(nb):
            ri, rj = divmod(r, m)
            for (pi, pj) in powers:
                c = ((ri + pi) % l) * m + ((rj + pj) % m)
                mat[r][c] ^= 1
        return mat

    A, B = circ(a_powers), circ(b_powers)
    # HX = [A | B], HZ = [B^T | A^T]
    HX = [A[r] + B[r] for r in range(nb)]
    HZ = [[B[c][r] for c in range(nb)] + [A[c][r] for c in range(nb)] for r in range(nb)]
    return HX, HZ


def write_mtx(path, rows):
    nr, nc = len(rows), len(rows[0])
    nnz = sum(sum(r) for r in rows)
    with open(path, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{nr} {nc} {nnz}\n")
        for i, r in enumerate(rows):
            for j, v in enumerate(r):
                if v:
                    fh.write(f"{i + 1} {j + 1} 1\n")


def write_qecc(path, HX, HZ, n, k, d):
    with open(path, "w") as fh:
        fh.write(f"{n},{k},{d}\n")
        for r in HX:                        # X-type stabilisers
            fh.write("".join(map(str, r)) + "0" * n + "\n")
        for r in HZ:                        # Z-type stabilisers
            fh.write("0" * n + "".join(map(str, r)) + "\n")
        fh.write("\n")


def rank2(rows):
    """GF(2) rank, to report the real k."""
    m = [r[:] for r in rows]
    piv, r = 0, 0
    ncols = len(m[0]) if m else 0
    for c in range(ncols):
        sel = next((i for i in range(r, len(m)) if m[i][c]), None)
        if sel is None:
            continue
        m[r], m[sel] = m[sel], m[r]
        for i in range(len(m)):
            if i != r and m[i][c]:
                m[i] = [a ^ b for a, b in zip(m[i], m[r])]
        r += 1
        piv += 1
    return piv


def main():
    kind = sys.argv[1]
    if kind == "tanner":
        HX, HZ, prefix, d = read_mtx(sys.argv[2]), read_mtx(sys.argv[3]), sys.argv[4], int(sys.argv[5])
        n = len(HX[0])
    else:
        l, m, prefix, d = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
        # the [[72,12,6]] family: A = x^3 + y + y^2, B = y^3 + x + x^2
        HX, HZ = bb_matrices(l, m, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)])
        n = 2 * l * m

    k = n - rank2(HX) - rank2(HZ)
    write_qecc(f"{prefix}.txt", HX, HZ, n, k, d)
    write_mtx(f"{prefix}_HX.mtx", HX)
    write_mtx(f"{prefix}_HZ.mtx", HZ)
    print(f"{prefix}: [[{n},{k},{d}]]  HX {len(HX)}x{n}  HZ {len(HZ)}x{n}")


if __name__ == "__main__":
    main()
