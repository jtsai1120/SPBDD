#include <algorithm>
// Scaling probe: run SPBDD's distance() over the codetables.de QECC benchmark
// from the codeDistance package, checking against the distance the file states.
//
//   file format:   n,k,d
//                  <2n binary chars>   [X block | Z block]
//                  ...
//                  <blank line>
#include <spbdd/spbdd.hpp>

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
    std::vector<std::string> stabilizers;   // as Pauli strings
    bool                     css = false;
};

// One row of the two-block matrix becomes one Pauli string.
std::string row_to_pauli(const std::string &bits, int n)
{
    std::string s(static_cast<std::size_t>(n), 'I');
    for (int q = 0; q < n; ++q) {
        const bool x = bits[static_cast<std::size_t>(q)] == '1';
        const bool z = bits[static_cast<std::size_t>(n + q)] == '1';
        s[static_cast<std::size_t>(q)] = pauli_to_char(make_pauli(x, z));
    }
    return s;
}

// A generator is CSS if it is all-X or all-Z; a code is CSS if all of them are.
bool row_is_css(const std::string &p)
{
    bool has_x = false, has_z = false;
    for (char c : p) {
        if (c == 'X') has_x = true;
        if (c == 'Z') has_z = true;
        if (c == 'Y') return false;
    }
    return !(has_x && has_z);
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

    for (Code &c : out) {
        c.css = true;
        for (const std::string &p : c.stabilizers)
            if (!row_is_css(p)) { c.css = false; break; }
    }
    return out;
}

} // namespace

int main(int argc, char **argv)
{
    const std::string path      = argc > 1 ? argv[1] : "QECC_sample.txt";
    const int         n_max     = argc > 2 ? std::atoi(argv[2]) : 16;
    const double      time_stop = argc > 3 ? std::atof(argv[3]) : 30.0;

    const std::vector<Code> codes = load(path);
    std::printf("loaded %zu codes from %s\n\n", codes.size(), path.c_str());
    std::printf("%-14s %-8s %6s %10s %12s %10s\n", "code", "type", "d", "spbdd d", "seconds",
                "nodes");

    for (const Code &c : codes) {
        if (c.n > n_max || c.stabilizers.empty()) continue;

        const auto t0 = Clock::now();
        int         got = -1;
        std::size_t nodes = 0;
        try {
            PauliSpace     sp(c.n);
            StabilizerCode code(sp, c.stabilizers);
            const PauliSet logical = code.normalizer() - code.group();
            nodes = logical.node_count();
            got   = logical.min_weight();
        } catch (const std::exception &e) {
            std::printf("[[%d,%d,%d]] threw: %s\n", c.n, c.k, c.d, e.what());
            continue;
        }
        const double secs = std::chrono::duration<double>(Clock::now() - t0).count();

        char name[32];
        std::snprintf(name, sizeof name, "[[%d,%d,%d]]", c.n, c.k, c.d);
        std::printf("%-14s %-8s %6d %10d %12.3f %10zu%s\n", name, c.css ? "CSS" : "non-CSS", c.d,
                    got, secs, nodes, got == c.d ? "" : "   <- MISMATCH");
        std::fflush(stdout);

        if (secs > time_stop) {
            std::printf("\nstopping: %.1f s exceeds the %.1f s budget at n=%d\n", secs, time_stop,
                        c.n);
            break;
        }
    }
    return 0;
}
