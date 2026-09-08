// Where does distance() spend its time at large d, and does advancing the
// weight bound incrementally help?
//
//   V1  what ships now: weight_at_most(w), rebuilt for each w
//   V2  weight_exactly(w) instead -- a thinner set to conjoin the checks into,
//       and equivalent while w is scanned upwards
//   V3  V2 with the checks conjoined smallest-support-first, to keep the
//       intermediate diagrams small
//
// V1 also reports where its time goes.
#include <spbdd/spbdd.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numeric>
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
        if (line.empty()) {
            if (found && !c.gens.empty()) break;
            continue;
        }
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

int weight_of(const std::string &p)
{
    return static_cast<int>(std::count_if(p.begin(), p.end(), [](char ch) { return ch != 'I'; }));
}

double secs(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// V1 -- what ships, with the phases timed.
int v1(PauliSpace &sp, const StabilizerCode &code, double &t_weight, double &t_checks,
       double &t_partners)
{
    t_weight = t_checks = t_partners = 0;
    for (int w = 1; w <= sp.n_qubits(); ++w) {
        auto      t0 = Clock::now();
        PauliSet  cand = sp.weight_at_most(w);
        t_weight += secs(t0);

        t0 = Clock::now();
        for (const std::string &g : code.stabilizers()) {
            cand &= sp.commuting_with(g);
            if (cand.is_empty()) break;
        }
        t_checks += secs(t0);
        if (cand.is_empty()) continue;

        t0 = Clock::now();
        for (const std::string &z : code.logical_z())
            if (!(cand & sp.anticommuting_with(z)).is_empty()) { t_partners += secs(t0); return w; }
        for (const std::string &x : code.logical_x())
            if (!(cand & sp.anticommuting_with(x)).is_empty()) { t_partners += secs(t0); return w; }
        t_partners += secs(t0);
    }
    return -1;
}

// V2 / V3 -- exactly-w, optionally with the checks ordered by support.
int v2(PauliSpace &sp, const StabilizerCode &code, bool order_by_support)
{
    std::vector<std::string> gens = code.stabilizers();
    if (order_by_support)
        std::stable_sort(gens.begin(), gens.end(),
                         [](const std::string &a, const std::string &b) {
                             return weight_of(a) < weight_of(b);
                         });

    for (int w = 1; w <= sp.n_qubits(); ++w) {
        PauliSet cand = sp.weight_exactly(w);
        for (const std::string &g : gens) {
            cand &= sp.commuting_with(g);
            if (cand.is_empty()) break;
        }
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
        if (c.gens.empty()) { std::printf("n=%s not found\n", argv[i]); continue; }

        std::printf("\n[[%d,%d,%d]]\n", c.n, c.k, c.d);

        double tw = 0, tc = 0, tp = 0, s1 = 0, s2 = 0, s3 = 0;
        int    d1 = 0, d2 = 0, d3 = 0;

        { PauliSpace sp(c.n); StabilizerCode co(sp, c.gens);
          auto t0 = Clock::now(); d1 = v1(sp, co, tw, tc, tp); s1 = secs(t0); }
        { PauliSpace sp(c.n); StabilizerCode co(sp, c.gens);
          auto t0 = Clock::now(); d2 = v2(sp, co, false); s2 = secs(t0); }
        { PauliSpace sp(c.n); StabilizerCode co(sp, c.gens);
          auto t0 = Clock::now(); d3 = v2(sp, co, true); s3 = secs(t0); }

        std::printf("  V1 at_most   %8.3f s   d=%d   [weight %5.1f%%  checks %5.1f%%  partners %5.1f%%]\n",
                    s1, d1, 100 * tw / s1, 100 * tc / s1, 100 * tp / s1);
        std::printf("  V2 exactly   %8.3f s   d=%d   %.2fx\n", s2, d2, s1 / s2);
        std::printf("  V3 +ordered  %8.3f s   d=%d   %.2fx\n", s3, d3, s1 / s3);
        std::fflush(stdout);
    }
    return 0;
}
