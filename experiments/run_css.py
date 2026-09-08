"""Same codes, but through the package's CSS path -- which splits the code into
its X and Z halves and solves two classical problems. That is what these tools
are built for, and it is the fair comparison on a CSS code like the toric one.

  usage: run_css.py <QECC file> <nMin> <nMax> <timeout s> <method>...
"""
import signal
import sys
import time

from codedistance import *                               # noqa: F401,F403
from codedistance.code_library import codeTables2Dict


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def load_codes(path):
    out, params, rows = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                if params and rows:
                    out.append((params, rows))
                params, rows = None, []
                continue
            if "," in line:
                if params and rows:
                    out.append((params, rows))
                    rows = []
                params = [int(x) for x in line.split(",")]
            elif params is not None:
                rows.append([int(c) for c in line])
    if params and rows:
        out.append((params, rows))
    return out


def main():
    path, n_min, n_max = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    budget, methods    = float(sys.argv[4]), sys.argv[5:]

    print(f"{'code':<14} {'d':>4} {'method':<16} {'got':>5} {'seconds':>10}  note", flush=True)
    signal.signal(signal.SIGALRM, _alarm)

    for (n, k, d), rows in load_codes(path):
        if n < n_min or n > n_max or not rows:
            continue
        S = codeTables2Dict([n, k, d], rows, "QECC")["S"]
        try:
            SX, SZ = CSSSplit(S)
        except Exception as exc:                          # noqa: BLE001
            print(f"[[{n},{k},{d}]] not CSS: {exc}")
            continue

        for method in methods:
            params = {"LOCheck": 0, "verbose": False, "nThreads": 1, "maxTime": budget}
            t0, note, got = time.perf_counter(), "", -1
            signal.setitimer(signal.ITIMER_REAL, budget)
            try:
                got = int(CSScodeDistance(SX, SZ, method=method, params=params)["d"])
            except Timeout:
                note = f"timeout >{budget:.0f}s"
            except Exception as exc:                      # noqa: BLE001
                note = f"{type(exc).__name__}: {str(exc)[:70]}"
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            secs = time.perf_counter() - t0
            if not note and got != d:
                note = f"got {got}, expected {d}"
            print(f"{f'[[{n},{k},{d}]]':<14} {d:>4} {method:<16} {got:>5} {secs:>10.3f}  {note}",
                  flush=True)


if __name__ == "__main__":
    main()
