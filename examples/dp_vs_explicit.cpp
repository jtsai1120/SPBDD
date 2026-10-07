// ===========================================================================
//  dp_vs_explicit: the DP of StabilizerCode::find_inequivalent_pair against an
//  explicit baseline, in time and in peak memory.
//
//      dp_vs_explicit <scenario...> <method>
//
//    scenarios (rotated surface code of distance d, k = 1)
//      ball   <d> <t>              every Pauli of weight <= t
//      region <d> <rows> <cols>    every Pauli supported on a rows x cols block (4^(rows*cols) of them)
//      regiont <d> <rows> <cols> <t>
//                                  every Pauli of weight <= t inside that block
//      cosets <d> <K> <G> <W>      union of K random cosets (G random generators of weight <= W)
//
//    methods
//      dp                          StabilizerCode::find_inequivalent_pair, PairMethod::Dp
//      explicit <cap_MB> [P0]      enumerate every element, keep a hash table sigma -> set of
//                                  lambda values seen, and stop when some sigma has two. The
//                                  table is capped at <cap_MB>; if it would exceed the cap the
//                                  enumeration is repeated with the syndrome space split into
//                                  P = 2, 4, 8, ... slices, one slice per pass
//                                  (starting from P0 when given).
//
//  The explicit baseline is told how the set is built (it enumerates by the recipe, one
//  XOR of precomputed vectors per element), which is the favourable case for it: the DP
//  receives only the diagram of the set. Both are single-threaded.
//  Output: verdict, wall time, peak resident memory of the whole process (the DP is run with the
//  BuDDy table of 2^22 nodes used by every other benchmark, which is part of its resident memory).
//
//  Build:  g++ -std=c++17 -O2 -Iinclude examples/dp_vs_explicit.cpp build/libspbdd.a -lm
// ===========================================================================

#include <spbdd/spbdd.hpp>

#include <sys/resource.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

using namespace spbdd;
using u64 = std::uint64_t;

static double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static double rss_mb()
{
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    return static_cast<double>(ru.ru_maxrss) / 1024.0;   // Linux: kilobytes
}

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

// ---------------------------------------------------------------------------
// Explicit baseline. pi = (sigma, lambda) packed in W words: sigma in bits 0..r-1,
// lambda (two bits, k = 1) in bits r, r+1. A table entry holds sigma in its key bits
// and, in bits r..r+3, the set of lambda values seen so far (bit r+lambda).
// ---------------------------------------------------------------------------
template <int W>
struct Pi {
    std::array<u64, W> w{};
    Pi &operator^=(const Pi &o)
    {
        for (int i = 0; i < W; ++i) w[i] ^= o.w[i];
        return *this;
    }
};

template <int W>
class Explicit {
public:
    Explicit(int r, std::size_t cap_bytes) : rw_(r / 64), rb_(r % 64), cap_bytes_(cap_bytes)
    {
        if (rb_ > 60) { std::fprintf(stderr, "bit layout unsupported for r=%d\n", r); std::exit(2); }
        lam_mask_ = (u64{3}) << rb_;
        mask_bits_ = (u64{15}) << rb_;
    }

    // One pass over the whole enumeration, handling only the sigma's of slice `pass` of `P`.
    // Returns: 0 = finished, no pair; 1 = pair found; 2 = table exceeded the cap.
    template <class Enumerate>
    int run_pass(const Enumerate &enumerate, unsigned P, unsigned pass)
    {
        P_ = P; pass_ = pass; found_ = false; over_ = false; count_ = 0; elements_ = 0;
        cap_ = 1u << 16;
        tab_.assign(cap_ * W, 0);
        enumerate([&](const Pi<W> &pi) { return visit(pi); });
        distinct_total_ += count_;
        peak_table_bytes_ = std::max(peak_table_bytes_, tab_.size() * sizeof(u64));
        if (over_) return 2;
        return found_ ? 1 : 0;
    }

    u64         distinct_total() const { return distinct_total_; }
    std::size_t peak_table_bytes() const { return peak_table_bytes_; }
    void        reset_stats() { distinct_total_ = 0; peak_table_bytes_ = 0; }

private:
    static u64 mix(u64 x)
    {
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
        return x;
    }

    bool visit(const Pi<W> &pi)   // returns false to stop the enumeration
    {
        const unsigned lam = static_cast<unsigned>((pi.w[rw_] >> rb_) & 3u);
        std::array<u64, W> key = pi.w;
        key[rw_] &= ~lam_mask_;
        u64 h = 0x9e3779b97f4a7c15ULL;
        for (int i = 0; i < W; ++i) h = mix(h ^ key[i]) + 0x9e3779b97f4a7c15ULL * (i + 1);
        if (P_ > 1 && static_cast<unsigned>((h >> 40) % P_) != pass_) return true;
        if (count_ * 10 >= cap_ * 7) {
            if (static_cast<std::size_t>(cap_) * 2 * W * sizeof(u64) > cap_bytes_) { over_ = true; return false; }
            grow();
        }
        std::size_t i = h & (cap_ - 1);
        for (;;) {
            u64 *e = &tab_[i * W];
            if ((e[rw_] & mask_bits_) == 0) {   // empty slot
                for (int j = 0; j < W; ++j) e[j] = key[j];
                e[rw_] |= u64{1} << (rb_ + lam);
                ++count_;
                return true;
            }
            bool same = true;
            for (int j = 0; j < W; ++j) {
                const u64 ej = (j == rw_) ? (e[j] & ~mask_bits_) : e[j];
                if (ej != key[j]) { same = false; break; }
            }
            if (same) {
                e[rw_] |= u64{1} << (rb_ + lam);
                const u64 m = (e[rw_] & mask_bits_) >> rb_;
                if (m & (m - 1)) { found_ = true; return false; }   // two lambda values for one sigma
                return true;
            }
            i = (i + 1) & (cap_ - 1);
        }
    }

    void grow()
    {
        std::vector<u64> old;
        old.swap(tab_);
        const std::size_t oc = cap_;
        cap_ *= 2;
        tab_.assign(cap_ * W, 0);
        for (std::size_t s = 0; s < oc; ++s) {
            const u64 *e = &old[s * W];
            if ((e[rw_] & mask_bits_) == 0) continue;
            std::array<u64, W> key;
            for (int j = 0; j < W; ++j) key[j] = (j == rw_) ? (e[j] & ~mask_bits_) : e[j];
            u64 h = 0x9e3779b97f4a7c15ULL;
            for (int i = 0; i < W; ++i) h = mix(h ^ key[i]) + 0x9e3779b97f4a7c15ULL * (i + 1);
            std::size_t i = h & (cap_ - 1);
            while ((tab_[i * W + rw_] & mask_bits_) != 0) i = (i + 1) & (cap_ - 1);
            for (int j = 0; j < W; ++j) tab_[i * W + j] = e[j];
        }
    }

    int         rw_, rb_;
    u64         lam_mask_, mask_bits_;
    std::size_t cap_bytes_;
    unsigned    P_ = 1, pass_ = 0;
    bool        found_ = false, over_ = false;
    std::size_t cap_ = 0;
    u64         count_ = 0, elements_ = 0, distinct_total_ = 0;
    std::size_t peak_table_bytes_ = 0;
    std::vector<u64> tab_;
};

// An enumeration is a function that calls `emit(pi)` once per element and stops when emit returns false.
template <int W>
using Emit = std::function<bool(const Pi<W> &)>;

template <int W>
static bool ball_rec(const std::vector<std::array<Pi<W>, 3>> &eff, int t, std::size_t start, int depth, const Pi<W> &cur,
                     const std::function<bool(const Pi<W> &)> &emit)
{
    if (!emit(cur)) return false;
    if (depth == t) return true;
    for (std::size_t q = start; q < eff.size(); ++q)
        for (int b = 0; b < 3; ++b) {
            Pi<W> nxt = cur;
            nxt ^= eff[q][static_cast<std::size_t>(b)];
            if (!ball_rec<W>(eff, t, q + 1, depth + 1, nxt, emit)) return false;
        }
    return true;
}

template <int W>
static bool gray(const Pi<W> &base, const std::vector<Pi<W>> &gens, const std::function<bool(const Pi<W> &)> &emit)
{
    Pi<W> cur = base;
    if (!emit(cur)) return false;
    const u64 total = u64{1} << gens.size();
    for (u64 i = 1; i < total; ++i) {
        cur ^= gens[static_cast<std::size_t>(__builtin_ctzll(i))];
        if (!emit(cur)) return false;
    }
    return true;
}

template <int W>
static int run_explicit(const StabilizerCode &code, int n, const std::string &scenario, const std::vector<std::string> &args,
                        std::size_t cap_mb, unsigned P0, double base_rss, double &secs_total, double &secs_final)
{
    const int r = code.n_stabilizers();
    auto pi_of = [&](const std::string &s) {
        Pi<W> p;
        const auto sy = code.syndrome(s);
        const auto sg = code.logical_signature(s);
        for (int i = 0; i < r; ++i) if (sy[static_cast<std::size_t>(i)]) p.w[static_cast<std::size_t>(i / 64)] |= u64{1} << (i % 64);
        for (std::size_t j = 0; j < sg.size(); ++j)
            if (sg[j]) p.w[static_cast<std::size_t>((r + static_cast<int>(j)) / 64)] |= u64{1} << ((r + static_cast<int>(j)) % 64);
        return p;
    };
    auto single = [&](int q, char c) { std::string s(static_cast<std::size_t>(n), 'I'); s[static_cast<std::size_t>(q)] = c; return pi_of(s); };

    // Build the recipe.
    std::function<bool(const Emit<W> &)> enumerate;
    std::vector<std::array<Pi<W>, 3>> eff;
    std::vector<Pi<W>>                gens;
    std::vector<std::pair<Pi<W>, std::vector<Pi<W>>>> cosets;
    int t = 0;
    const int d = std::atoi(args[0].c_str());
    if (scenario == "ball" || scenario == "regiont") {
        std::vector<int> qs;
        if (scenario == "ball") {
            t = std::atoi(args[1].c_str());
            for (int q = 0; q < n; ++q) qs.push_back(q);
        } else {
            const int rows = std::atoi(args[1].c_str()), cols = std::atoi(args[2].c_str());
            t = std::atoi(args[3].c_str());
            for (int a = 0; a < rows; ++a) for (int b = 0; b < cols; ++b) qs.push_back(a * d + b);
        }
        for (int q : qs) {
            const Pi<W> x = single(q, 'X'), z = single(q, 'Z');
            Pi<W> y = x; y ^= z;
            eff.push_back({x, z, y});
        }
        enumerate = [&](const Emit<W> &emit) { return ball_rec<W>(eff, t, 0, 0, Pi<W>{}, emit); };
    } else if (scenario == "region") {
        const int rows = std::atoi(args[1].c_str()), cols = std::atoi(args[2].c_str());
        for (int a = 0; a < rows; ++a)
            for (int b = 0; b < cols; ++b) { gens.push_back(single(a * d + b, 'X')); gens.push_back(single(a * d + b, 'Z')); }
        enumerate = [&](const Emit<W> &emit) { return gray<W>(Pi<W>{}, gens, emit); };
    } else if (scenario == "cosets") {
        const int K = std::atoi(args[1].c_str()), G = std::atoi(args[2].c_str()), Wt = std::atoi(args[3].c_str());
        unsigned long long state = 88172645463325252ULL;   // same generator as pair_benchmark
        auto rnd = [&]() { state ^= state << 13; state ^= state >> 7; state ^= state << 17; return state; };
        auto pauli = [&](int max_weight) {
            std::string s(static_cast<std::size_t>(n), 'I');
            const int   w = 1 + static_cast<int>(rnd() % static_cast<unsigned>(max_weight));
            for (int i = 0; i < w; ++i) s[rnd() % static_cast<unsigned>(n)] = "XYZ"[rnd() % 3];
            return s;
        };
        for (int c = 0; c < K; ++c) {
            std::vector<std::string> gs;
            for (int i = 0; i < G; ++i) gs.push_back(pauli(Wt));
            const std::string base = pauli(Wt);
            std::vector<Pi<W>> gp;
            for (auto &g : gs) gp.push_back(pi_of(g));
            cosets.emplace_back(pi_of(base), gp);
        }
        enumerate = [&](const Emit<W> &emit) {
            for (auto &c : cosets) if (!gray<W>(c.first, c.second, emit)) return false;
            return true;
        };
    } else {
        std::fprintf(stderr, "unknown scenario\n");
        return -1;
    }

    // Passes, doubling P until the table fits under the cap.
    int      verdict = 0;
    double   total0  = now();
    unsigned P       = P0;
    std::printf("  explicit: cap %zu MB\n", cap_mb);
    for (;;) {
        Explicit<W> ex(r, cap_mb << 20);
        const double a = now();
        bool aborted = false;
        verdict = 0;
        for (unsigned pass = 0; pass < P; ++pass) {
            const int rc = ex.run_pass([&](const auto &emit) { return enumerate(emit); }, P, pass);
            if (rc == 2) { aborted = true; break; }
            if (rc == 1) { verdict = 1; break; }
        }
        const double b = now();
        if (aborted) {
            std::printf("    P=%-3u table exceeded the cap after %.1f s; doubling P\n", P, b - a);
            P *= 2;
            continue;
        }
        secs_final = b - a;
        std::printf("    P=%-3u finished in %.1f s; distinct sigma seen %.4g, largest table %.0f MB\n", P, secs_final,
                    static_cast<double>(ex.distinct_total()), static_cast<double>(ex.peak_table_bytes()) / (1 << 20));
        break;
    }
    secs_total = now() - total0;
    (void)base_rss;
    return verdict;
}

int main(int argc, char **argv)
{
    auto usage = [&]() {
        std::fprintf(stderr,
                     "usage: %s ball <d> <t>                 (dp | explicit <cap_MB>)\n"
                     "       %s region <d> <rows> <cols>     (dp | explicit <cap_MB>)\n"
                     "       %s regiont <d> <rows> <cols> <t> (dp | explicit <cap_MB>)\n"
                     "       %s cosets <d> <K> <G> <W>       (dp | explicit <cap_MB>)\n",
                     argv[0], argv[0], argv[0], argv[0]);
        return 2;
    };
    if (argc < 4) return usage();
    const std::string scenario = argv[1];
    const int         nargs = scenario == "ball" ? 2 : scenario == "region" ? 3 : scenario == "regiont" ? 4 : scenario == "cosets" ? 4 : -1;
    if (nargs < 0 || argc < 2 + nargs + 1) return usage();
    std::vector<std::string> args;
    for (int i = 0; i < nargs; ++i) args.push_back(argv[2 + i]);
    const std::string method = argv[2 + nargs];
    const int         d      = std::atoi(args[0].c_str());
    const int         n      = d * d;

    // The DP gets the table size used for every other benchmark. The explicit baseline
    // uses the library only to compute syndromes, so it gets a small one and its
    // resident memory is then essentially its own hash table.
    ManagerConfig cfg;
    const bool    is_dp = (method == "dp");
    cfg.initial_nodes = is_dp ? (1u << 22) : (1u << 14);
    cfg.cache_size    = is_dp ? (1u << 20) : (1u << 12);
    PauliSpace           sp(n, cfg);
    const StabilizerCode code(sp, surface(d));
    const double         base = rss_mb();
    std::printf("%s", scenario.c_str());
    for (auto &a : args) std::printf(" %s", a.c_str());
    std::printf("   [%s]   start-up resident memory %.0f MB\n", method.c_str(), base);

    if (method == "dp") {
        PauliSet errors = sp.empty();
        if (scenario == "ball") {
            errors = sp.weight_at_most(std::atoi(args[1].c_str()));
        } else if (scenario == "region" || scenario == "regiont") {
            const int rows = std::atoi(args[1].c_str()), cols = std::atoi(args[2].c_str());
            std::vector<int> qs;
            for (int a = 0; a < rows; ++a) for (int b = 0; b < cols; ++b) qs.push_back(a * d + b);
            errors = sp.supported_on(qs);   // identity outside the block
            if (scenario == "regiont") errors = errors & sp.weight_at_most(std::atoi(args[3].c_str()), qs);
        } else {
            const int K = std::atoi(args[1].c_str()), G = std::atoi(args[2].c_str()), Wt = std::atoi(args[3].c_str());
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
                for (int i = 0; i < G; ++i) gens.push_back(pauli(Wt));
                errors |= sp.coset_of(pauli(Wt), gens);
            }
        }
        std::printf("  set: %.4g elements, %zu nodes\n", errors.size(), errors.node_count());
        const double t0 = now();
        const auto   f  = code.find_inequivalent_pair(errors, StabilizerCode::PairMethod::Dp);
        const double dt = now() - t0;
        std::printf("RESULT dp        pair=%s  time %.3f s  peak memory %.0f MB  (BuDDy table %zu nodes)\n",
                    f ? "yes" : "no", dt, rss_mb(), sp.manager().peak_nodes());
        return 0;
    }
    if (method == "explicit" && argc >= 4 + nargs) {
        const std::size_t cap_mb = static_cast<std::size_t>(std::atoll(argv[3 + nargs]));
        const unsigned    P0     = argc > 4 + nargs ? static_cast<unsigned>(std::atoi(argv[4 + nargs])) : 1u;
        const int r = code.n_stabilizers();
        double secs_total = 0, secs_final = 0;
        int verdict = 0;
        if (r + 4 <= 128)      verdict = run_explicit<2>(code, n, scenario, args, cap_mb, P0, base, secs_total, secs_final);
        else if (r + 4 <= 192) verdict = run_explicit<3>(code, n, scenario, args, cap_mb, P0, base, secs_total, secs_final);
        else                   verdict = run_explicit<4>(code, n, scenario, args, cap_mb, P0, base, secs_total, secs_final);
        std::printf("RESULT explicit  pair=%s  time %.3f s (final configuration), %.3f s (including aborted attempts)  peak memory %.0f MB\n",
                    verdict == 1 ? "yes" : "no", secs_final, secs_total, rss_mb());
        return 0;
    }
    return usage();
}
