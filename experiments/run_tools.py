"""Run the codeDistance non-CSS-capable methods over the codetables.de QECC
benchmark, on this machine, so the comparison against SPBDD is like for like.

Each method fills in its own parameter defaults through the package's
setDefaultParams, so the settings match the ones the paper ran with.

  usage: run_tools.py <QECC file> <nMin> <nMax> <per-code timeout s> <method>...
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
    """codetables.de QECC file -> [((n, k, d), rows)] with rows in [X|Z] form."""
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
    path    = sys.argv[1]
    n_min   = int(sys.argv[2])
    n_max   = int(sys.argv[3])
    budget  = float(sys.argv[4])
    methods = sys.argv[5:]

    codes = load_codes(path)
    print(f"{'code':<14} {'type':<8} {'d':>4} {'method':<20} {'got':>5} {'seconds':>10}  note",
          flush=True)
    signal.signal(signal.SIGALRM, _alarm)

    for (n, k, d), rows in codes:
        if n < n_min or n > n_max or not rows:
            continue
        S = codeTables2Dict([n, k, d], rows, "QECC")["S"]

        # non-CSS iff some generator mixes X and Z on one qubit
        css = all(not (any(r[q] for q in range(n)) and any(r[n + q] for q in range(n)))
                  for r in rows)

        for method in methods:
            # Only the keys codeDistance itself reads; the rest are filled in
            # by each method's own setDefaultParams call.
            params = {"LOCheck": 0, "verbose": False, "nThreads": 1, "maxTime": budget}
            t0, note, got = time.perf_counter(), "", -1
            signal.setitimer(signal.ITIMER_REAL, budget)
            try:
                got = int(codeDistance(S, None, 2, method=method, params=params)["d"])
            except Timeout:
                note = f"timeout >{budget:.0f}s"
            except Exception as exc:                      # noqa: BLE001
                note = f"{type(exc).__name__}: {str(exc)[:70]}"
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            secs = time.perf_counter() - t0

            if not note and got != d:
                note = "MISMATCH (over)" if got > d else "upper bound only"
            print(f"{f'[[{n},{k},{d}]]':<14} {'CSS' if css else 'non-CSS':<8} {d:>4} "
                  f"{method:<20} {got:>5} {secs:>10.3f}  {note}", flush=True)


if __name__ == "__main__":
    main()
