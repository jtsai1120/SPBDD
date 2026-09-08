// Does constraining early help?
//
// (a) is what distance() does today: build N(S), subtract S, then bisect on the
//     weight bound. The whole logical set is materialised first.
//
// (b) is the Brouwer-Zimmermann idea in decision-diagram form: start from the
//     weight bound, which is small while w is small, and conjoin the parity
//     checks one at a time. Nothing ever has to represent N(S) in full.
//
// Both compute the same number; only the order of the conjunction differs.
#include <spbdd/spbdd.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace spbdd;
using Clock = std::chrono::steady_clock;

namespace {

struct Code {
    int                      n, k, d;
    std::vector<std::string> stabilizers;
};

std::string row_to_pauli(const std::string &bits, int n)
{
    std::string s(static_cast<std::size_t>(n), 'I');
    for (int q = 0; q < n; ++q)
        s[static_cast<std::size_t>(q)] =
            pauli_to_char(make_pauli(bits[static_cast<std::size_t>(q)] == '1',
                                     bits[static_cast<std::size_t>(n + q)] == '1'));
    return s;
}

std::vector<Code> load(const std::string &path)
{
    std::ifstream     in(path);
    std::vector<Code> out;
    std::string       line;
    Code              cur;
    bool              open = false;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
        if (line.empty()) {
            if (open) { out.push_back(cur); open = false; }
            continue;
        }
        if (line.find(',') != std::string::npos) {
            if (open) out.push_back(cur);
            cur = Code{};
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream ss(line);
            ss >> cur.n >> cur.k >> cur.d;
            open = true;
            continue;
        }
        if (open && static_cast<int>(line.size()) == 2 * cur.n)
            cur.stabilizers.push_back(row_to_pauli(line, cur.n));
    }
    if (open) out.push_back(cur);
    return out;
}

// (b) the weight bound comes first and the checks are conjoined into it.
int distance_constrained(PauliSpace &sp, const StabilizerCode &code, int n, std::size_t &peak)
{
    // "the logical signature is not zero" = anticommutes with at least one
    // partner. Built once; it is a union of 2k parity checks.
    PauliSet nonzero_signature = sp.empty();
    for (const std::string &z : code.logical_z()) nonzero_signature |= sp.anticommuting_with(z);
    for (const std::string &x : code.logical_x()) nonzero_signature |= sp.anticommuting_with(x);

    for (int w = 1; w <= n; ++w) {
        PauliSet acc = sp.weight_at_most(w);
        for (const std::string &g : code.stabilizers()) {
            acc &= sp.commuting_with(g);
            if (acc.is_empty()) break;          // nothing of this weight survives
        }
        if (acc.is_empty()) continue;
        acc &= nonzero_signature;
        peak = std::max(peak, acc.node_count());
        if (!acc.is_empty()) return w;
    }
    return -1;
}

} // namespace

int main(int argc, char **argv)
{
    const std::string path  = argv[1];
    const int         n_min = std::atoi(argv[2]);
    const int         n_max = std::atoi(argv[3]);

    std::printf("%-13s %4s | %10s %11s | %10s %11s | %s\n", "code", "d", "(a) secs", "(a) nodes",
                "(b) secs", "(b) nodes", "speedup");
    for (const Code &c : load(path)) {
        if (c.n < n_min || c.n > n_max || c.stabilizers.empty()) continue;

        double      sa = 0, sb = 0;
        int         da = -1, db = -1;
        std::size_t na = 0, nb = 0;

        {   // (a) materialise N(S) - S, then bisect
            PauliSpace     sp(c.n);
            StabilizerCode code(sp, c.stabilizers);
            const auto     t0 = Clock::now();
            const PauliSet logical = code.normalizer() - code.group();
            na = logical.node_count();
            da = logical.min_weight();
            sa = std::chrono::duration<double>(Clock::now() - t0).count();
        }
        {   // (b) constrain first
            PauliSpace     sp(c.n);
            StabilizerCode code(sp, c.stabilizers);
            const auto     t0 = Clock::now();
            db = distance_constrained(sp, code, c.n, nb);
            sb = std::chrono::duration<double>(Clock::now() - t0).count();
        }

        char name[32];
        std::snprintf(name, sizeof name, "[[%d,%d,%d]]", c.n, c.k, c.d);
        std::printf("%-13s %4d | %10.3f %11zu | %10.3f %11zu | %6.1fx%s\n", name, c.d, sa, na, sb,
                    nb, sb > 0 ? sa / sb : 0.0,
                    (da == c.d && db == c.d) ? "" : "   <- WRONG");
        std::fflush(stdout);
    }
    return 0;
}
