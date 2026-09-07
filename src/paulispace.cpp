// ===========================================================================
//  PauliSpace -- construction of Pauli operator sets.
//
//  Three of the builders here would be embarrassing if written from the
//  definition, and the way they avoid that is the point of the library:
//
//    * generated_by does not multiply generators together. It computes the
//      null space of the generator matrix and turns each basis vector into one
//      parity check, so the cost is independent of the 2^rank group size.
//
//    * weight_* does not union C(n,w) cubes. It carries a counter along the
//      qubits, which is O(n*w).
//
//    * commuting_with is a single parity check, not a filter over candidates.
// ===========================================================================

#include "spbdd/paulispace.hpp"
#include "spbdd/pauliset.hpp"

#include <stdexcept>
#include <string>

namespace spbdd {
namespace {

using Row = std::vector<std::uint8_t>;   // a vector over GF(2), length 2n

void require_qubit(const PauliSpace &sp, int q, const char *who)
{
    if (q < 0 || q >= sp.n_qubits())
        throw std::out_of_range(std::string(who) + ": qubit " + std::to_string(q) +
                                " outside a space of " + std::to_string(sp.n_qubits()));
}

// Basis of { h : g . h == 0 for every row g of M }, i.e. the orthogonal
// complement of the row space. Plain Gaussian elimination to reduced row
// echelon form, then one basis vector per free column.
std::vector<Row> null_space(std::vector<Row> m, int n_cols)
{
    std::vector<int> pivot_of_row;
    int              row = 0;

    for (int col = 0; col < n_cols && row < static_cast<int>(m.size()); ++col) {
        int found = -1;
        for (int r = row; r < static_cast<int>(m.size()); ++r)
            if (m[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)]) { found = r; break; }
        if (found < 0) continue;

        std::swap(m[static_cast<std::size_t>(row)], m[static_cast<std::size_t>(found)]);
        for (int r = 0; r < static_cast<int>(m.size()); ++r) {
            if (r == row) continue;
            if (!m[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)]) continue;
            for (int c = 0; c < n_cols; ++c)
                m[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] ^=
                    m[static_cast<std::size_t>(row)][static_cast<std::size_t>(c)];
        }
        pivot_of_row.push_back(col);
        ++row;
    }

    std::vector<char> is_pivot(static_cast<std::size_t>(n_cols), 0);
    for (int c : pivot_of_row) is_pivot[static_cast<std::size_t>(c)] = 1;

    // Setting one free variable to 1 and the others to 0 forces each pivot
    // variable to the corresponding coefficient in its row.
    std::vector<Row> basis;
    for (int free = 0; free < n_cols; ++free) {
        if (is_pivot[static_cast<std::size_t>(free)]) continue;
        Row h(static_cast<std::size_t>(n_cols), 0);
        h[static_cast<std::size_t>(free)] = 1;
        for (std::size_t i = 0; i < pivot_of_row.size(); ++i)
            h[static_cast<std::size_t>(pivot_of_row[i])] = m[i][static_cast<std::size_t>(free)];
        basis.push_back(std::move(h));
    }
    return basis;
}

// A Pauli string as a GF(2) vector under the xvar/zvar convention.
Row to_row(const std::vector<Pauli> &p)
{
    Row v(2 * p.size(), 0);
    for (std::size_t q = 0; q < p.size(); ++q) {
        v[2 * q]     = has_x(p[q]) ? 1 : 0;
        v[2 * q + 1] = has_z(p[q]) ? 1 : 0;
    }
    return v;
}

std::uint8_t dot(const Row &a, const Row &b)
{
    std::uint8_t s = 0;
    for (std::size_t i = 0; i < a.size(); ++i) s ^= static_cast<std::uint8_t>(a[i] & b[i]);
    return s;
}

std::vector<int> support(const Row &h)
{
    std::vector<int> vars;
    for (std::size_t i = 0; i < h.size(); ++i)
        if (h[i]) vars.push_back(static_cast<int>(i));
    return vars;
}

std::vector<int> all_qubits(int n)
{
    std::vector<int> qs(static_cast<std::size_t>(n));
    for (int q = 0; q < n; ++q) qs[static_cast<std::size_t>(q)] = q;
    return qs;
}

// The symplectic form against a fixed operator is a linear functional, so
// "commutes with p" is one parity check: the x variable of a qubit appears
// exactly when p has a z there, and vice versa.
std::vector<int> commutation_vars(const std::vector<Pauli> &p)
{
    std::vector<int> vars;
    for (std::size_t q = 0; q < p.size(); ++q) {
        if (has_z(p[q])) vars.push_back(static_cast<int>(2 * q));       // xvar(q)
        if (has_x(p[q])) vars.push_back(static_cast<int>(2 * q + 1));   // zvar(q)
    }
    return vars;
}

} // namespace

// ===========================================================================
//  Construction and geometry
// ===========================================================================

PauliSpace::PauliSpace(int n_qubits, const ManagerConfig &cfg)
{
    if (n_qubits < 0) throw std::invalid_argument("PauliSpace: negative qubit count");
    st_ = std::make_shared<State>(n_qubits, cfg);
}

int PauliSpace::n_qubits() const { return st_->n; }
int PauliSpace::n_vars() const { return 2 * st_->n; }

Manager &PauliSpace::manager() const { return st_->mgr; }

void PauliSpace::grow_to(int n_qubits)
{
    if (n_qubits < st_->n)
        throw std::invalid_argument("grow_to: a PauliSpace can only grow");
    if (n_qubits == st_->n) return;
    st_->mgr.ensure_var_count(2 * n_qubits);
    st_->n = n_qubits;
}

// ===========================================================================
//  Primitives
// ===========================================================================

PauliSet PauliSpace::wrap(Bdd f) const
{
    if (f.valid() && f.manager() != &st_->mgr)
        throw std::logic_error("wrap: the Bdd belongs to a different Manager");
    return PauliSet(*this, std::move(f));
}

PauliSet PauliSpace::empty() const { return wrap(st_->mgr.constant(false)); }

// Every one of the 4^n assignments is a valid Pauli string, so the universe is
// simply the constant true and the complement of a set needs no correction.
PauliSet PauliSpace::all() const { return wrap(st_->mgr.constant(true)); }

PauliSet PauliSpace::identity() const
{
    Bdd f = st_->mgr.constant(true);
    for (int v = 0; v < n_vars(); ++v) f &= st_->mgr.literal(v, false);
    return wrap(std::move(f));
}

PauliSet PauliSpace::literal(int var, bool positive) const
{
    if (var < 0 || var >= n_vars())
        throw std::out_of_range("literal: variable outside the space");
    return wrap(st_->mgr.literal(var, positive));
}

PauliSet PauliSpace::cube(const std::vector<std::pair<int, bool>> &assignment) const
{
    Bdd f = st_->mgr.constant(true);
    for (const auto &a : assignment) {
        if (a.first < 0 || a.first >= n_vars())
            throw std::out_of_range("cube: variable outside the space");
        f &= st_->mgr.literal(a.first, a.second);
    }
    return wrap(std::move(f));
}

PauliSet PauliSpace::parity(const std::vector<int> &vars, bool target) const
{
    // [ x1 ^ ... ^ xk == target ] is the negation of the XOR against target,
    // which is the same as folding the literals into the constant !target.
    Bdd f = st_->mgr.constant(!target);
    for (int v : vars) {
        if (v < 0 || v >= n_vars())
            throw std::out_of_range("parity: variable outside the space");
        f ^= st_->mgr.literal(v, true);
    }
    return wrap(std::move(f));
}

// ===========================================================================
//  Explicitly listed elements
// ===========================================================================

PauliSet PauliSpace::from(const std::vector<Pauli> &p) const
{
    if (static_cast<int>(p.size()) != st_->n)
        throw std::invalid_argument("from: operator has the wrong number of qubits");

    Bdd f = st_->mgr.constant(true);
    for (int q = 0; q < st_->n; ++q) {
        f &= st_->mgr.literal(xvar(q), has_x(p[static_cast<std::size_t>(q)]));
        f &= st_->mgr.literal(zvar(q), has_z(p[static_cast<std::size_t>(q)]));
    }
    return wrap(std::move(f));
}

PauliSet PauliSpace::from(const std::string &s) const
{
    return from(parse_pauli_string(s, st_->n, "from"));
}

PauliSet PauliSpace::from_list(const std::vector<std::string> &list) const
{
    Bdd f = st_->mgr.constant(false);
    for (const std::string &s : list) f |= from(s).bdd();
    return wrap(std::move(f));
}

// ===========================================================================
//  Subgroups
// ===========================================================================
//
//  Let W be the GF(2) span of the generators. Multiplying them out would mean
//  2^rank elements. Instead, use that the form is non-degenerate, so
//  (W^perp)^perp = W:
//
//      v is in W   <=>   h . v == 0 for every h in a basis of W^perp
//
//  W^perp is the null space of the generator matrix, which one Gaussian
//  elimination produces. Each basis vector then becomes a single parity check
//  of about 2n nodes, and the cost never sees the size of the group.

PauliSet PauliSpace::generated_by(const std::vector<std::string> &generators) const
{
    std::vector<Row> rows;
    rows.reserve(generators.size());
    for (const std::string &g : generators)
        rows.push_back(to_row(parse_pauli_string(g, st_->n, "generated_by")));

    Bdd f = st_->mgr.constant(true);
    for (const Row &h : null_space(rows, n_vars())) f &= parity(support(h), false).bdd();
    return wrap(std::move(f));
}

PauliSet PauliSpace::coset_of(const std::string &base,
                              const std::vector<std::string> &generators) const
{
    const Row base_row = to_row(parse_pauli_string(base, st_->n, "coset_of"));

    std::vector<Row> rows;
    rows.reserve(generators.size());
    for (const std::string &g : generators)
        rows.push_back(to_row(parse_pauli_string(g, st_->n, "coset_of")));

    // Same checks as the subgroup, but each one now has to reproduce the value
    // the base element gives it.
    Bdd f = st_->mgr.constant(true);
    for (const Row &h : null_space(rows, n_vars()))
        f &= parity(support(h), dot(h, base_row) != 0).bdd();
    return wrap(std::move(f));
}

// ===========================================================================
//  Support
// ===========================================================================

PauliSet PauliSpace::supported_on(const std::vector<int> &qubits) const
{
    std::vector<char> allowed(static_cast<std::size_t>(st_->n), 0);
    for (int q : qubits) {
        require_qubit(*this, q, "supported_on");
        allowed[static_cast<std::size_t>(q)] = 1;
    }

    Bdd f = st_->mgr.constant(true);
    for (int q = 0; q < st_->n; ++q) {
        if (allowed[static_cast<std::size_t>(q)]) continue;
        f &= st_->mgr.literal(xvar(q), false);
        f &= st_->mgr.literal(zvar(q), false);
    }
    return wrap(std::move(f));
}

PauliSet PauliSpace::identity_on(const std::vector<int> &qubits) const
{
    Bdd f = st_->mgr.constant(true);
    for (int q : qubits) {
        require_qubit(*this, q, "identity_on");
        f &= st_->mgr.literal(xvar(q), false);
        f &= st_->mgr.literal(zvar(q), false);
    }
    return wrap(std::move(f));
}

PauliSet PauliSpace::pauli_at(int qubit, Pauli p) const
{
    require_qubit(*this, qubit, "pauli_at");
    Bdd f = st_->mgr.literal(xvar(qubit), has_x(p));
    f &= st_->mgr.literal(zvar(qubit), has_z(p));
    return wrap(std::move(f));
}

PauliSet PauliSpace::pauli_at(int qubit, char p) const
{
    return pauli_at(qubit, pauli_from_char(p));
}

PauliSet PauliSpace::non_identity_at(int qubit) const
{
    require_qubit(*this, qubit, "non_identity_at");
    Bdd f = st_->mgr.literal(xvar(qubit), true) | st_->mgr.literal(zvar(qubit), true);
    return wrap(std::move(f));
}

PauliSet PauliSpace::matching(const std::string &pattern) const
{
    if (static_cast<int>(pattern.size()) != st_->n)
        throw std::invalid_argument("matching: pattern has the wrong number of qubits");

    Bdd f = st_->mgr.constant(true);
    for (int q = 0; q < st_->n; ++q) {
        const char c = pattern[static_cast<std::size_t>(q)];
        if (c == '*' || c == '?' || c == '.') continue;   // this qubit is free
        const Pauli p = pauli_from_char(c);
        f &= st_->mgr.literal(xvar(q), has_x(p));
        f &= st_->mgr.literal(zvar(q), has_z(p));
    }
    return wrap(std::move(f));
}

// ===========================================================================
//  Weight
// ===========================================================================
//
//  Taken from the definition this is a union of C(n,w) cubes. Instead keep, as
//  a small array of diagrams, "exactly c of the qubits seen so far are not the
//  identity", and extend it one qubit at a time:
//
//      next[c] = (cur[c] & identity_q) | (cur[c-1] & non_identity_q)
//
//  Counts above the upper bound are dropped as they appear, so the work is
//  O(n*hi) operations on diagrams of O(n*hi) nodes.

PauliSet PauliSpace::weight_between(int lo, int hi, const std::vector<int> &qubits) const
{
    for (int q : qubits) require_qubit(*this, q, "weight_between");

    const int n_counted = static_cast<int>(qubits.size());
    if (lo < 0) lo = 0;
    if (hi > n_counted) hi = n_counted;
    if (lo > hi) return empty();

    std::vector<Bdd> cur;
    cur.reserve(static_cast<std::size_t>(hi) + 1);
    cur.push_back(st_->mgr.constant(true));
    for (int c = 1; c <= hi; ++c) cur.push_back(st_->mgr.constant(false));

    for (int q : qubits) {
        const Bdd x           = st_->mgr.literal(xvar(q), true);
        const Bdd z           = st_->mgr.literal(zvar(q), true);
        const Bdd is_identity = (!x) & (!z);
        const Bdd is_other    = x | z;

        std::vector<Bdd> next;
        next.reserve(static_cast<std::size_t>(hi) + 1);
        for (int c = 0; c <= hi; ++c) {
            Bdd stay = cur[static_cast<std::size_t>(c)] & is_identity;
            if (c > 0) stay |= cur[static_cast<std::size_t>(c - 1)] & is_other;
            next.push_back(std::move(stay));
        }
        cur = std::move(next);
    }

    Bdd f = st_->mgr.constant(false);
    for (int c = lo; c <= hi; ++c) f |= cur[static_cast<std::size_t>(c)];
    return wrap(std::move(f));
}

PauliSet PauliSpace::weight_between(int lo, int hi) const
{
    return weight_between(lo, hi, all_qubits(st_->n));
}

PauliSet PauliSpace::weight_exactly(int w) const { return weight_between(w, w); }
PauliSet PauliSpace::weight_at_most(int w) const { return weight_between(0, w); }

PauliSet PauliSpace::weight_exactly(int w, const std::vector<int> &qubits) const
{
    return weight_between(w, w, qubits);
}

PauliSet PauliSpace::weight_at_most(int w, const std::vector<int> &qubits) const
{
    return weight_between(0, w, qubits);
}

// ===========================================================================
//  Commutation
// ===========================================================================

PauliSet PauliSpace::commuting_with(const std::string &p) const
{
    return parity(commutation_vars(parse_pauli_string(p, st_->n, "commuting_with")), false);
}

PauliSet PauliSpace::anticommuting_with(const std::string &p) const
{
    return parity(commutation_vars(parse_pauli_string(p, st_->n, "anticommuting_with")), true);
}

PauliSet PauliSpace::commuting_with_all(const std::vector<std::string> &ps) const
{
    Bdd f = st_->mgr.constant(true);
    for (const std::string &p : ps) f &= commuting_with(p).bdd();
    return wrap(std::move(f));
}

} // namespace spbdd
