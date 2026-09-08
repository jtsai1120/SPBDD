"""Emit toric codes in the codetables.de QECC format, so SPBDD and the
codeDistance tools read byte-identical input.

Toric code on a d x d torus: qubits on edges, n = 2d^2, k = 2, distance d.
Every check touches four edges that are neighbours on the lattice, which is the
sparse/local case -- the opposite of the codetables.de codes.

  usage: gen_toric.py <d_min> <d_max> > toric.txt
"""
import sys


def emit(d):
    n = 2 * d * d

    def h(i, j):                      # horizontal edge (i,j)
        return (i % d) * d + (j % d)

    def v(i, j):                      # vertical edge (i,j)
        return d * d + (i % d) * d + (j % d)

    rows = []
    # vertex operators carry X; the last one is dependent on the rest
    for i in range(d):
        for j in range(d):
            if i == d - 1 and j == d - 1:
                continue
            x = [0] * n
            for e in (h(i, j), h(i, j - 1), v(i, j), v(i - 1, j)):
                x[e] ^= 1
            rows.append(x + [0] * n)

    # plaquette operators carry Z; likewise one is dependent
    for i in range(d):
        for j in range(d):
            if i == d - 1 and j == d - 1:
                continue
            z = [0] * n
            for e in (h(i, j), h(i + 1, j), v(i, j), v(i, j + 1)):
                z[e] ^= 1
            rows.append([0] * n + z)

    print(f"{n},2,{d}")
    for r in rows:
        print("".join(str(b) for b in r))
    print()


for d in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    emit(d)
