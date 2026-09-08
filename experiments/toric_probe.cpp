// Is SPBDD's ceiling n, or is it structure?
//
// The codetables.de codes are dense and unstructured -- the worst case for a
// decision diagram, whose size depends on how far apart in the variable order a
// check reaches. The toric code is the opposite: every check touches four
// edges that are adjacent in a row-major indexing. If the ceiling is bandwidth
// rather than n, this should go much further.
//
//   toric code on a d x d torus: n = 2d^2 qubits (edges), k = 2, distance d
#include <spbdd/spbdd.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace spbdd;
using Clock = std::chrono::steady_clock;

namespace {

// Horizontal edge (i,j) -> i*d + j ; vertical edge (i,j) -> d*d + i*d + j
int h_edge(int d, int i, int j) { return ((i % d + d) % d) * d + ((j % d + d) % d); }
int v_edge(int d, int i, int j) { return d * d + ((i % d + d) % d) * d + ((j % d + d) % d); }

std::string put(int n, const std::vector<int> &qubits, char letter)
{
    std::string s(static_cast<std::size_t>(n), 'I');
    for (int q : qubits) s[static_cast<std::size_t>(q)] = letter;
    return s;
}

// Vertex operators carry X, plaquette operators carry Z. One of each family is
// dependent on the rest, so the last one of each is dropped.
std::vector<std::string> toric(int d)
{
    const int                n = 2 * d * d;
    std::vector<std::string> out;

    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j) {
            if (i == d - 1 && j == d - 1) continue;   // dependent
            out.push_back(put(n, {h_edge(d, i, j), h_edge(d, i, j - 1),
                                  v_edge(d, i, j), v_edge(d, i - 1, j)}, 'X'));
        }

    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j) {
            if (i == d - 1 && j == d - 1) continue;   // dependent
            out.push_back(put(n, {h_edge(d, i, j), h_edge(d, i + 1, j),
                                  v_edge(d, i, j), v_edge(d, i, j + 1)}, 'Z'));
        }
    return out;
}

} // namespace

int main(int argc, char **argv)
{
    const int    d_max     = argc > 1 ? std::atoi(argv[1]) : 6;
    const double time_stop = argc > 2 ? std::atof(argv[2]) : 120.0;

    std::printf("%-16s %5s %5s %8s %10s %12s %12s\n", "toric code", "n", "k", "want d", "spbdd d",
                "seconds", "nodes");

    for (int d = 2; d <= d_max; ++d) {
        const int  n    = 2 * d * d;
        const auto gens = toric(d);

        const auto  t0    = Clock::now();
        int         got   = -1, k = -1;
        std::size_t nodes = 0;
        try {
            PauliSpace     sp(n);
            StabilizerCode code(sp, gens);
            k = code.n_logical();
            const PauliSet logical = code.normalizer() - code.group();
            nodes = logical.node_count();
            got   = logical.min_weight();
        } catch (const std::exception &e) {
            std::printf("d=%d n=%d threw: %s\n", d, n, e.what());
            break;
        }
        const double secs = std::chrono::duration<double>(Clock::now() - t0).count();

        char name[32];
        std::snprintf(name, sizeof name, "d=%d", d);
        std::printf("%-16s %5d %5d %8d %10d %12.3f %12zu%s\n", name, n, k, d, got, secs, nodes,
                    got == d ? "" : "   <- MISMATCH");
        std::fflush(stdout);

        if (secs > time_stop) {
            std::printf("\nstopping: %.1f s exceeds the %.1f s budget at n=%d\n", secs, time_stop, n);
            break;
        }
    }
    return 0;
}
