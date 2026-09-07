#pragma once

// ===========================================================================
//  A brute-force model of the library, for the tests to compare against.
//
//  Every function here is written straight from the definition, with sets of
//  strings and explicit loops -- exactly the implementations the library goes
//  out of its way to avoid. That is the point: a test that compares a clever
//  implementation against a second clever implementation only proves the two
//  agree. Comparing against the definition proves the definition is met.
//
//  The cost is 4^n, so callers keep n at 3 or 4. `span` is 2^k in the number
//  of generators.
// ===========================================================================

#include <spbdd/spbdd.hpp>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace spbdd_ref {

using Strings = std::set<std::string>;

// All 4^n Pauli strings on n qubits.
inline Strings universe(int n)
{
    Strings out{""};
    for (int q = 0; q < n; ++q) {
        Strings next;
        for (const std::string &s : out)
            for (char c : std::string("IXZY")) next.insert(s + c);
        out = std::move(next);
    }
    return out;
}

// --- set algebra, done the obvious way -------------------------------------

inline Strings unite(const Strings &a, const Strings &b)
{
    Strings r = a;
    r.insert(b.begin(), b.end());
    return r;
}

inline Strings intersect(const Strings &a, const Strings &b)
{
    Strings r;
    for (const std::string &s : a)
        if (b.count(s)) r.insert(s);
    return r;
}

inline Strings subtract(const Strings &a, const Strings &b)
{
    Strings r;
    for (const std::string &s : a)
        if (!b.count(s)) r.insert(s);
    return r;
}

inline Strings sym_diff(const Strings &a, const Strings &b)
{
    return unite(subtract(a, b), subtract(b, a));
}

inline Strings complement(const Strings &a, int n) { return subtract(universe(n), a); }

// --- the builders ----------------------------------------------------------

// The span of the generators, formed by actually multiplying out all 2^k
// products.
inline Strings span(const std::vector<std::string> &generators, int n)
{
    Strings   out;
    const int k = static_cast<int>(generators.size());
    for (int mask = 0; mask < (1 << k); ++mask) {
        std::vector<spbdd::Pauli> acc(static_cast<std::size_t>(n), spbdd::Pauli::I);
        for (int i = 0; i < k; ++i)
            if (mask & (1 << i))
                acc = spbdd::pauli_mul(
                    acc, spbdd::parse_pauli_string(generators[static_cast<std::size_t>(i)]));
        out.insert(spbdd::pauli_string_to_text(acc));
    }
    return out;
}

inline Strings coset(const std::string &base, const std::vector<std::string> &generators, int n)
{
    Strings out;
    for (const std::string &s : span(generators, n))
        out.insert(spbdd::pauli_string_to_text(
            spbdd::pauli_mul(spbdd::parse_pauli_string(base), spbdd::parse_pauli_string(s))));
    return out;
}

// Filter the whole universe by weight, counting only the listed qubits (empty
// means all of them).
inline Strings weight_between(int n, int lo, int hi, const std::vector<int> &qubits = {})
{
    Strings out;
    for (const std::string &s : universe(n)) {
        const auto p = spbdd::parse_pauli_string(s);
        int        w = 0;
        if (qubits.empty()) {
            w = spbdd::pauli_weight(p);
        } else {
            for (int q : qubits)
                if (p[static_cast<std::size_t>(q)] != spbdd::Pauli::I) ++w;
        }
        if (w >= lo && w <= hi) out.insert(s);
    }
    return out;
}

inline Strings commuting_with(const std::string &fixed, int n)
{
    Strings    out;
    const auto b = spbdd::parse_pauli_string(fixed);
    for (const std::string &s : universe(n))
        if (spbdd::pauli_commute(spbdd::parse_pauli_string(s), b)) out.insert(s);
    return out;
}

// Multiply every element by every Pauli supported on the listed qubits.
inline Strings forget(const Strings &b, const std::vector<int> &qubits, int n)
{
    Strings out;
    for (const std::string &s : b) {
        const auto base = spbdd::parse_pauli_string(s);
        for (const std::string &g : universe(n)) {
            const auto gp = spbdd::parse_pauli_string(g);
            bool       ok = true;
            for (int q = 0; q < n && ok; ++q) {
                const bool listed = std::find(qubits.begin(), qubits.end(), q) != qubits.end();
                if (!listed && gp[static_cast<std::size_t>(q)] != spbdd::Pauli::I) ok = false;
            }
            if (ok) out.insert(spbdd::pauli_string_to_text(spbdd::pauli_mul(base, gp)));
        }
    }
    return out;
}

// Rewrite every element to be the identity on the listed qubits.
inline Strings reset(const Strings &b, const std::vector<int> &qubits)
{
    Strings out;
    for (const std::string &s : b) {
        auto p = spbdd::parse_pauli_string(s);
        for (int q : qubits) p[static_cast<std::size_t>(q)] = spbdd::Pauli::I;
        out.insert(spbdd::pauli_string_to_text(p));
    }
    return out;
}

// --- conjugation -----------------------------------------------------------
//
// A Clifford conjugation is defined by where it sends the generators X_q and
// Z_q. So rather than repeat the library's coordinate formulas, decompose each
// operator into its X and Z generators and multiply the stated images back
// together -- a different route to the same map, which is the point.
//
// `images` lists the image of X, then Z, for each qubit the gate touches, as a
// local Pauli string on those qubits. For CX on {c,t} that is
// {"XX", "ZI", "IX", "ZZ"}: X_c -> X_c X_t, Z_c -> Z_c, X_t -> X_t,
// Z_t -> Z_c Z_t.

inline std::vector<spbdd::Pauli> embed(const std::string &local,
                                       const std::vector<int> &qubits, int n)
{
    std::vector<spbdd::Pauli> out(static_cast<std::size_t>(n), spbdd::Pauli::I);
    const auto                p = spbdd::parse_pauli_string(local);
    for (std::size_t i = 0; i < qubits.size(); ++i)
        out[static_cast<std::size_t>(qubits[i])] = p[i];
    return out;
}

inline Strings conjugate(const Strings &b, const std::vector<int> &qubits,
                         const std::vector<std::string> &images, int n)
{
    Strings out;
    for (const std::string &s : b) {
        const auto p = spbdd::parse_pauli_string(s);

        // Qubits the gate does not touch keep their generators, so copy them.
        std::vector<spbdd::Pauli> r(static_cast<std::size_t>(n), spbdd::Pauli::I);
        for (int q = 0; q < n; ++q)
            if (std::find(qubits.begin(), qubits.end(), q) == qubits.end())
                r[static_cast<std::size_t>(q)] = p[static_cast<std::size_t>(q)];

        // The touched ones are rebuilt from the images of their generators.
        for (std::size_t i = 0; i < qubits.size(); ++i) {
            const spbdd::Pauli here = p[static_cast<std::size_t>(qubits[i])];
            if (spbdd::has_x(here)) r = spbdd::pauli_mul(r, embed(images[2 * i], qubits, n));
            if (spbdd::has_z(here)) r = spbdd::pauli_mul(r, embed(images[2 * i + 1], qubits, n));
        }
        out.insert(spbdd::pauli_string_to_text(r));
    }
    return out;
}

// --- crossing back over ----------------------------------------------------

// The same set as a PauliSet, built element by element, so a test can compare
// it with the library's answer using PauliSet::operator==.
inline spbdd::PauliSet as_set(const spbdd::PauliSpace &sp, const Strings &s)
{
    return sp.from_list(std::vector<std::string>(s.begin(), s.end()));
}

inline std::vector<std::string> as_vector(const Strings &s)
{
    return std::vector<std::string>(s.begin(), s.end());
}

} // namespace spbdd_ref
