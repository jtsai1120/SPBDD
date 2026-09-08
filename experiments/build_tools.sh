#!/bin/bash
# Build the reference tools this comparison runs against, without root.
#
#   codeDistance   the package from arXiv:2603.22532 (BZDistMW, connectedClusterMW,
#                  pySATDist, QDistEvol)
#   m4riCC/m4riRW  dist-m4ri, which the paper recommends for quantum CSS codes
#
# The distributions want `sudo apt install python3-pip libm4ri-dev`; neither is
# needed. pip is bootstrapped into ~/.local and m4ri is built into a prefix
# under /tmp.
set -eu

echo "=== pip, into ~/.local ==="
if ! python3 -m pip --version >/dev/null 2>&1; then
    curl -sSLo /tmp/get-pip.py https://bootstrap.pypa.io/get-pip.py
    python3 /tmp/get-pip.py --user -q
fi
python3 -m pip --version

echo "=== codeDistance ==="
python3 -m pip install --user -q codedistance python-sat
python3 -c "import codedistance; print('codedistance ok')"
# the package's own examples/ needs Python 3.12 f-strings, so the harnesses
# here talk to codedistance directly rather than importing runTest.py
git clone -q --depth 1 https://github.com/m-webster/codeDistancePYPI.git /tmp/codeDistance 2>/dev/null || true

echo "=== m4ri (release tarball: it ships configure, the git tree does not) ==="
if [ ! -f /tmp/m4ri-install/lib/libm4ri.a ]; then
    cd /tmp
    curl -sSLO https://github.com/malb/m4ri/releases/download/20260122/m4ri-20260122.tar.gz
    tar xzf m4ri-20260122.tar.gz
    cd m4ri-20260122
    ./configure --prefix=/tmp/m4ri-install --disable-shared >/dev/null
    # same maintainer-mode trap as CUDD: git-fresh timestamps make make try to
    # regenerate with an automake nobody has
    make -j"$(nproc)" ACLOCAL=: AUTOCONF=: AUTOMAKE=: AUTOHEADER=: >/dev/null
    make install >/dev/null
fi
ls -la /tmp/m4ri-install/lib/libm4ri.a

echo "=== dist-m4ri ==="
[ -d /tmp/dist-m4ri ] || git clone -q --depth 1 https://github.com/QEC-pages/dist-m4ri.git /tmp/dist-m4ri
cd /tmp/dist-m4ri/src
# the makefile does not thread LDFLAGS into the link line, so -L rides on CFLAGS
make all CFLAGS="-O3 -march=native -I/tmp/m4ri-install/include -L/tmp/m4ri-install/lib" >/dev/null
./dist_m4ri 2>&1 | head -1

echo
echo "all tools ready"
