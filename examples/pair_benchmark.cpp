// ===========================================================================
//  pair_benchmark: time StabilizerCode::find_inequivalent_pair, DP vs compose.
//
//      pair_benchmark ball   <d> <t> <dp|compose>
//          rotated surface code of distance d, error set = every Pauli of
//          weight <= t.
//      pair_benchmark cosets <d> <K> <G> <W> <dp|compose>
//          the same code, error set = union of K random cosets, each of G
//          random generators of weight <= W (a set with no structure to lean
//          on). Deterministic: fixed seed.
//
//  Build:  g++ -std=c++17 -O2 -Iinclude examples/pair_benchmark.cpp build/libspbdd.a -lm
// ===========================================================================

#include <spbdd/spbdd.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace spbdd;

static double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Rotated surface code on a d x d grid: weight-4 plaquettes in the bulk,
// weight-2 ones on alternate boundary edges. n = d^2, n - 1 generators, k = 1.
static std::vector<std::string> surface(int d)
{
    const int                n = d * d;
    std::vector<std::string> g;
    for (int i = 0; i <= d; ++i)
        for (int j = 0; j <= d; ++j) {
            std::vector<int> sup;
            for (int dr = -1; dr <= 0; ++dr)
                for (int dc = -1; dc <= 0; ++dc) {
                    const int r = i + dr, c = j + dc;
                    if (r >= 0 && r < d && c >= 0 && c < d) sup.push_back(r * d + c);
                }
            const bool x = ((i + j) % 2 == 0);
            if (sup.size() == 2) {
                const bool horizontal = (i == 0 || i == d);
                if (horizontal && !x) continue;
                if (!horizontal && x) continue;
            } else if (sup.size() != 4) {
                continue;
            }
            std::string s(static_cast<std::size_t>(n), 'I');
            for (int v : sup) s[static_cast<std::size_t>(v)] = x ? 'X' : 'Z';
            g.push_back(s);
        }
    return g;
}

int main(int argc, char **argv)
{
    auto usage = [&]() {
        std::fprintf(stderr,
                     "usage: %s ball <d> <t> <dp|compose>\n"
                     "       %s cosets <d> <K> <G> <W> <dp|compose>\n", argv[0], argv[0]);
        return 2;
    };
    if (argc < 5) return usage();

    const std::string mode = argv[1];
    const int         d    = std::atoi(argv[2]);
    const int         n    = d * d;
    const std::string method_name = (mode == "ball") ? argv[4] : (argc > 6 ? argv[6] : "dp");
    const auto        method = (method_name == "compose") ? StabilizerCode::PairMethod::Compose
                                                           : StabilizerCode::PairMethod::Dp;

    // BuDDy's node table is fixed at start-up and grows only when it fills, and
    // a table that is too small makes every method garbage-collect constantly,
    // so the numbers are only comparable between methods at one table size.
    ManagerConfig cfg;
    cfg.initial_nodes = 1u << 22;
    cfg.cache_size    = 1u << 20;
    PauliSpace           sp(n, cfg);
    const StabilizerCode code(sp, surface(d));

    PauliSet errors = sp.empty();
    if (mode == "ball") {
        errors = sp.weight_at_most(std::atoi(argv[3]));
    } else if (mode == "cosets" && argc >= 7) {
        const int K = std::atoi(argv[3]), G = std::atoi(argv[4]), W = std::atoi(argv[5]);
        unsigned long long state = 88172645463325252ULL;
        auto rnd = [&]() { state ^= state << 13; state ^= state >> 7; state ^= state << 17; return state; };
        auto pauli = [&](int max_weight) {
            std::string s(static_cast<std::size_t>(n), 'I');
            const int   w = 1 + static_cast<int>(rnd() % static_cast<unsigned>(max_weight));
            for (int i = 0; i < w; ++i) s[rnd() % static_cast<unsigned>(n)] = "XYZ"[rnd() % 3];
            return s;
        };
        for (int c = 0; c < K; ++c) {
            std::vector<std::string> gens;
            for (int i = 0; i < G; ++i) gens.push_back(pauli(W));
            errors |= sp.coset_of(pauli(W), gens);
        }
    } else {
        return usage();
    }

    std::printf("surface d=%d (n=%d)  errors: %.4g elements, %zu nodes\n", d, n, errors.size(), errors.node_count());
    const double t0    = now();
    const auto   found = code.find_inequivalent_pair(errors, method);
    std::printf("%-8s  inequivalent pair: %s   %.3f s\n", method_name == "compose" ? "compose" : "dp",
                found ? "yes" : "no", now() - t0);
    return 0;
}
