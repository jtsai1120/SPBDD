// The conjunction of the stabilizer checks is 99.9% of the time, so this looks
// at that phase only.
//
//   per-w breakdown  -- is the last w the whole cost, or is the scan the cost?
//   V4 sifting on    -- the intermediates are large, so let the package reorder
//   V5 interleaved   -- conjoin a check as soon as its variables are all seen
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
    std::vector<std::string> gens;
};

Code load(const std::string &path, int want)
{
    std::ifstream in(path);
    std::string   line;
    Code          c{};
    bool          found = false;
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
        if (line.empty()) { if (found && !c.gens.empty()) break; continue; }
        if (line.find(',') != std::string::npos) {
            if (found && !c.gens.empty()) break;
            std::string t = line;
            std::replace(t.begin(), t.end(), ',', ' ');
            std::istringstream ss(t);
            ss >> c.n >> c.k >> c.d;
            found = (c.n == want);
            c.gens.clear();
            continue;
        }
        if (found && static_cast<int>(line.size()) == 2 * c.n) {
            std::string p(static_cast<std::size_t>(c.n), 'I');
            for (int q = 0; q < c.n; ++q)
                p[static_cast<std::size_t>(q)] =
                    pauli_to_char(make_pauli(line[static_cast<std::size_t>(q)] == '1',
                                             line[static_cast<std::size_t>(c.n + q)] == '1'));
            c.gens.push_back(p);
        }
    }
    return c;
}

double secs(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

int run(PauliSpace &sp, const StabilizerCode &code, bool report_per_w, std::size_t *peak = nullptr)
{
    for (int w = 1; w <= sp.n_qubits(); ++w) {
        const auto t0   = Clock::now();
        PauliSet   cand = sp.weight_at_most(w);
        std::size_t biggest = cand.node_count();
        for (const std::string &g : code.stabilizers()) {
            cand &= sp.commuting_with(g);
            biggest = std::max(biggest, cand.node_count());
            if (cand.is_empty()) break;
        }
        const double t = secs(t0);
        if (peak) *peak = std::max(*peak, biggest);
        if (report_per_w)
            std::printf("      w=%-3d %8.3f s   biggest intermediate %9zu nodes\n", w, t, biggest);
        if (cand.is_empty()) continue;
        for (const std::string &z : code.logical_z())
            if (!(cand & sp.anticommuting_with(z)).is_empty()) return w;
        for (const std::string &x : code.logical_x())
            if (!(cand & sp.anticommuting_with(x)).is_empty()) return w;
    }
    return -1;
}

} // namespace

int main(int argc, char **argv)
{
    const std::string path = argv[1];
    for (int i = 2; i < argc; ++i) {
        const Code c = load(path, std::atoi(argv[i]));
        if (c.gens.empty()) continue;
        std::printf("\n[[%d,%d,%d]]  r=%zu checks\n", c.n, c.k, c.d, c.gens.size());

        double      base = 0, sift = 0;
        int         d0 = 0, d1 = 0;
        std::size_t p0 = 0, p1 = 0;

        {   std::printf("  per-w, reordering off:\n");
            PauliSpace sp(c.n); StabilizerCode co(sp, c.gens);
            auto t0 = Clock::now(); d0 = run(sp, co, true, &p0); base = secs(t0); }

        {   PauliSpace sp(c.n); StabilizerCode co(sp, c.gens);
            sp.manager().set_dynamic_reordering(true);
            auto t0 = Clock::now(); d1 = run(sp, co, false, &p1); sift = secs(t0); }

        std::printf("  total off  %8.3f s  d=%d  peak intermediate %9zu\n", base, d0, p0);
        std::printf("  sifting on %8.3f s  d=%d  peak intermediate %9zu   %.2fx\n", sift, d1, p1,
                    base / sift);
        std::fflush(stdout);
    }
    return 0;
}
