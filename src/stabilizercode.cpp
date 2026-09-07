// ===========================================================================
//  StabilizerCode -- the symplectic basis, the change of coordinates it
//  defines, and the inequivalent-pair decision procedure.
//
//  All of the linear algebra runs once, in the constructor. A query is then a
//  fixed number of diagram operations plus O(k), whatever the error set holds.
// ===========================================================================

#include "spbdd/stabilizercode.hpp"

#include <stdexcept>
#include <string>

namespace spbdd {
namespace {

using Row = std::vector<std::uint8_t>;

Row to_row(const std::vector<Pauli> &p)
{
    Row v(2 * p.size(), 0);
    for (std::size_t q = 0; q < p.size(); ++q) {
        v[2 * q]     = has_x(p[q]) ? 1 : 0;
        v[2 * q + 1] = has_z(p[q]) ? 1 : 0;
    }
    return v;
}

std::string row_to_string(const Row &v)
{
    std::vector<Pauli> p(v.size() / 2, Pauli::I);
    for (std::size_t q = 0; q < p.size(); ++q) p[q] = make_pauli(v[2 * q] != 0, v[2 * q + 1] != 0);
    return pauli_string_to_text(p);
}

// The symplectic form pairs the x part of one operator with the z part of the
// other, so it is the ordinary dot product against the partner with x and z
// exchanged. Every "read a coefficient by pairing with the partner" below is
// this function.
Row swap_xz(const Row &v)
{
    Row out(v.size(), 0);
    for (std::size_t q = 0; q * 2 + 1 < v.size(); ++q) {
        out[2 * q]     = v[2 * q + 1];
        out[2 * q + 1] = v[2 * q];
    }
    return out;
}

std::uint8_t dot(const Row &a, const Row &b)
{
    std::uint8_t s = 0;
    for (std::size_t i = 0; i < a.size(); ++i) s ^= static_cast<std::uint8_t>(a[i] & b[i]);
    return s;
}

std::uint8_t symplectic(const Row &a, const Row &b) { return dot(a, swap_xz(b)); }

void add_into(Row &target, const Row &source)
{
    for (std::size_t i = 0; i < target.size(); ++i) target[i] ^= source[i];
}

bool is_zero(const Row &v)
{
    for (std::uint8_t x : v)
        if (x) return false;
    return true;
}

// Reduce to row echelon form, dropping zero rows. The surviving rows are
// independent and span the same space, which is all the caller needs.
std::vector<Row> independent_rows(std::vector<Row> m, int n_cols)
{
    std::vector<Row> out;
    int              row = 0;
    for (int col = 0; col < n_cols && row < static_cast<int>(m.size()); ++col) {
        int found = -1;
        for (int rr = row; rr < static_cast<int>(m.size()); ++rr)
            if (m[static_cast<std::size_t>(rr)][static_cast<std::size_t>(col)]) { found = rr; break; }
        if (found < 0) continue;
        std::swap(m[static_cast<std::size_t>(row)], m[static_cast<std::size_t>(found)]);
        for (int rr = 0; rr < static_cast<int>(m.size()); ++rr)
            if (rr != row && m[static_cast<std::size_t>(rr)][static_cast<std::size_t>(col)])
                add_into(m[static_cast<std::size_t>(rr)], m[static_cast<std::size_t>(row)]);
        ++row;
    }
    for (int i = 0; i < row; ++i) out.push_back(m[static_cast<std::size_t>(i)]);
    return out;
}

// Reduced row echelon form of m, together with the transform e such that
// r = e * m, and the pivot column of each surviving row.
struct Reduction {
    std::vector<Row> r, e;
    std::vector<int> pivot;
};

Reduction reduce(std::vector<Row> m, int n_cols)
{
    const std::size_t rows = m.size();
    std::vector<Row>  e(rows, Row(rows, 0));
    for (std::size_t i = 0; i < rows; ++i) e[i][i] = 1;

    Reduction out;
    int       row = 0;
    for (int col = 0; col < n_cols && row < static_cast<int>(rows); ++col) {
        int found = -1;
        for (int rr = row; rr < static_cast<int>(rows); ++rr)
            if (m[static_cast<std::size_t>(rr)][static_cast<std::size_t>(col)]) { found = rr; break; }
        if (found < 0) continue;

        std::swap(m[static_cast<std::size_t>(row)], m[static_cast<std::size_t>(found)]);
        std::swap(e[static_cast<std::size_t>(row)], e[static_cast<std::size_t>(found)]);
        for (int rr = 0; rr < static_cast<int>(rows); ++rr) {
            if (rr == row) continue;
            if (!m[static_cast<std::size_t>(rr)][static_cast<std::size_t>(col)]) continue;
            add_into(m[static_cast<std::size_t>(rr)], m[static_cast<std::size_t>(row)]);
            add_into(e[static_cast<std::size_t>(rr)], e[static_cast<std::size_t>(row)]);
        }
        out.pivot.push_back(col);
        ++row;
    }
    for (int i = 0; i < row; ++i) {
        out.r.push_back(m[static_cast<std::size_t>(i)]);
        out.e.push_back(e[static_cast<std::size_t>(i)]);
    }
    return out;
}

// Basis of { v : row . v == 0 for every row of m }.
std::vector<Row> null_space(const std::vector<Row> &m, int n_cols)
{
    const Reduction  red = reduce(m, n_cols);
    std::vector<char> is_pivot(static_cast<std::size_t>(n_cols), 0);
    for (int c : red.pivot) is_pivot[static_cast<std::size_t>(c)] = 1;

    std::vector<Row> basis;
    for (int free = 0; free < n_cols; ++free) {
        if (is_pivot[static_cast<std::size_t>(free)]) continue;
        Row h(static_cast<std::size_t>(n_cols), 0);
        h[static_cast<std::size_t>(free)] = 1;
        for (std::size_t i = 0; i < red.pivot.size(); ++i)
            h[static_cast<std::size_t>(red.pivot[i])] = red.r[i][static_cast<std::size_t>(free)];
        basis.push_back(std::move(h));
    }
    return basis;
}

// Gauss-Jordan on [ m | I ]. m must be square and non-singular.
std::vector<Row> invert(const std::vector<Row> &m)
{
    const int        n = static_cast<int>(m.size());
    std::vector<Row> a = m;
    std::vector<Row> b(static_cast<std::size_t>(n), Row(static_cast<std::size_t>(n), 0));
    for (int i = 0; i < n; ++i) b[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1;

    for (int col = 0; col < n; ++col) {
        int found = -1;
        for (int rr = col; rr < n; ++rr)
            if (a[static_cast<std::size_t>(rr)][static_cast<std::size_t>(col)]) { found = rr; break; }
        if (found < 0) throw std::logic_error("StabilizerCode: the coordinate transform is singular");

        std::swap(a[static_cast<std::size_t>(col)], a[static_cast<std::size_t>(found)]);
        std::swap(b[static_cast<std::size_t>(col)], b[static_cast<std::size_t>(found)]);
        for (int rr = 0; rr < n; ++rr) {
            if (rr == col) continue;
            if (!a[static_cast<std::size_t>(rr)][static_cast<std::size_t>(col)]) continue;
            add_into(a[static_cast<std::size_t>(rr)], a[static_cast<std::size_t>(col)]);
            add_into(b[static_cast<std::size_t>(rr)], b[static_cast<std::size_t>(col)]);
        }
    }
    return b;
}

} // namespace

// ===========================================================================
//  Construction: the symplectic basis
// ===========================================================================

StabilizerCode::StabilizerCode(PauliSpace space, const std::vector<std::string> &generators)
    : sp_(std::move(space))
{
    const int n      = sp_.n_qubits();
    const int n_vars = sp_.n_vars();

    // --- the stabilizers ---------------------------------------------------
    std::vector<Row> given;
    given.reserve(generators.size());
    for (const std::string &s : generators)
        given.push_back(to_row(parse_pauli_string(s, n, "StabilizerCode")));

    for (std::size_t i = 0; i < given.size(); ++i)
        for (std::size_t j = i + 1; j < given.size(); ++j)
            if (symplectic(given[i], given[j]))
                throw std::invalid_argument("StabilizerCode: generators " + std::to_string(i) +
                                            " and " + std::to_string(j) + " do not commute");

    grow_ = independent_rows(given, n_vars);
    r_    = static_cast<int>(grow_.size());
    k_    = n - r_;

    // --- the destabilizers -------------------------------------------------
    // Wanted: d_i with <g_j, d_i> = delta_ji. Since <g_j, v> = swap_xz(g_j) . v
    // that is one linear system with the same matrix for every i, so reduce it
    // once and read off a particular solution per right-hand side.
    std::vector<Row> pairing;
    pairing.reserve(grow_.size());
    for (const Row &g : grow_) pairing.push_back(swap_xz(g));

    const Reduction red = reduce(pairing, n_vars);
    if (static_cast<int>(red.pivot.size()) != r_)
        throw std::logic_error("StabilizerCode: the stabilizer pairing map is not surjective");

    std::vector<Row> draft;
    for (int i = 0; i < r_; ++i) {
        // With every non-pivot coordinate zero, the reduced system forces the
        // pivot coordinates to be the transformed right-hand side.
        Row v(static_cast<std::size_t>(n_vars), 0);
        for (int j = 0; j < r_; ++j)
            v[static_cast<std::size_t>(red.pivot[static_cast<std::size_t>(j)])] =
                red.e[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
        draft.push_back(std::move(v));
    }

    // The drafts pair correctly with the stabilizers but not yet with each
    // other. Adding stabilizers fixes that without disturbing the pairing:
    // <d_i + sum c_ik g_k, d_j + ...> = <d_i,d_j> + c_ij + c_ji, so putting the
    // offending product into the upper triangle cancels it.
    drow_ = draft;
    for (int i = 0; i < r_; ++i)
        for (int j = i + 1; j < r_; ++j)
            if (symplectic(draft[static_cast<std::size_t>(i)], draft[static_cast<std::size_t>(j)]))
                add_into(drow_[static_cast<std::size_t>(i)], grow_[static_cast<std::size_t>(j)]);

    // --- the logical operators ---------------------------------------------
    // What commutes with every stabilizer and every destabilizer: a 2k
    // dimensional space on which the form is still non-degenerate.
    std::vector<Row> constraints;
    for (const Row &g : grow_) constraints.push_back(swap_xz(g));
    for (const Row &d : drow_) constraints.push_back(swap_xz(d));

    std::vector<Row> l = null_space(constraints, n_vars);

    // Symplectic Gram-Schmidt: pull out one hyperbolic pair at a time and make
    // everything left orthogonal to both, so later pairs cannot interfere.
    while (!l.empty()) {
        const Row u = l.front();
        l.erase(l.begin());
        if (is_zero(u)) continue;

        std::size_t partner = l.size();
        for (std::size_t i = 0; i < l.size(); ++i)
            if (symplectic(u, l[i])) { partner = i; break; }
        if (partner == l.size())
            throw std::logic_error("StabilizerCode: the logical space is degenerate");

        const Row w = l[partner];
        l.erase(l.begin() + static_cast<std::ptrdiff_t>(partner));

        for (Row &v : l) {
            if (symplectic(v, w)) add_into(v, u);
            if (symplectic(v, u)) add_into(v, w);
        }

        lxrow_.push_back(u);
        lzrow_.push_back(w);
    }

    if (static_cast<int>(lxrow_.size()) != k_)
        throw std::logic_error("StabilizerCode: found " + std::to_string(lxrow_.size()) +
                               " logical pairs, expected " + std::to_string(k_));

    for (const Row &v : grow_) g_.push_back(row_to_string(v));
    for (const Row &v : drow_) d_.push_back(row_to_string(v));
    for (const Row &v : lxrow_) lx_.push_back(row_to_string(v));
    for (const Row &v : lzrow_) lz_.push_back(row_to_string(v));

    // --- the coordinate transform ------------------------------------------
    // Row v of T reads new coordinate v out of an operator, and each is read by
    // pairing with the *partner* of the basis vector it is the coefficient of.
    std::vector<Row> t(static_cast<std::size_t>(n_vars), Row(static_cast<std::size_t>(n_vars), 0));
    for (int i = 0; i < r_; ++i) {
        t[static_cast<std::size_t>(a_var(i))] = swap_xz(drow_[static_cast<std::size_t>(i)]);
        t[static_cast<std::size_t>(b_var(i))] = swap_xz(grow_[static_cast<std::size_t>(i)]);
    }
    for (int j = 0; j < k_; ++j) {
        t[static_cast<std::size_t>(c_var(j))] = swap_xz(lzrow_[static_cast<std::size_t>(j)]);
        t[static_cast<std::size_t>(f_var(j))] = swap_xz(lxrow_[static_cast<std::size_t>(j)]);
    }
    tinv_ = invert(t);
}

// ===========================================================================
//  The sets
// ===========================================================================

PauliSet StabilizerCode::group() const { return sp_.generated_by(g_); }
PauliSet StabilizerCode::normalizer() const { return sp_.commuting_with_all(g_); }

int StabilizerCode::distance() const { return (normalizer() - group()).min_weight(); }

// ===========================================================================
//  Reading one operator
// ===========================================================================

std::vector<bool> StabilizerCode::syndrome(const std::string &e) const
{
    const Row          v = to_row(parse_pauli_string(e, sp_.n_qubits(), "syndrome"));
    std::vector<bool>  out;
    for (const Row &g : grow_) out.push_back(symplectic(v, g) != 0);
    return out;
}

std::vector<bool> StabilizerCode::logical_signature(const std::string &e) const
{
    const Row         v = to_row(parse_pauli_string(e, sp_.n_qubits(), "logical_signature"));
    std::vector<bool> out;
    for (const Row &z : lzrow_) out.push_back(symplectic(v, z) != 0);   // c
    for (const Row &x : lxrow_) out.push_back(symplectic(v, x) != 0);   // f
    return out;
}

// ===========================================================================
//  Reading a whole class
// ===========================================================================
//
//  "Anticommutes with g_i exactly when bit i is set" is a conjunction of parity
//  checks, so these are built in the original coordinates and never touch the
//  transform. That makes them an independent route to the same partition.

namespace {

PauliSet class_of(const PauliSpace &sp, const std::vector<std::string> &partners,
                  const std::vector<bool> &bits, const char *who)
{
    if (bits.size() != partners.size())
        throw std::invalid_argument(std::string(who) + ": expected " +
                                    std::to_string(partners.size()) + " bits, got " +
                                    std::to_string(bits.size()));

    PauliSet out = sp.all();
    for (std::size_t i = 0; i < partners.size(); ++i)
        out = out & (bits[i] ? sp.anticommuting_with(partners[i]) : sp.commuting_with(partners[i]));
    return out;
}

} // namespace

PauliSet StabilizerCode::with_syndrome(const std::vector<bool> &syndrome) const
{
    return class_of(sp_, g_, syndrome, "with_syndrome");
}

PauliSet StabilizerCode::with_logical_signature(const std::vector<bool> &signature) const
{
    // c is read against Zbar and f against Xbar, so the partners come in that
    // order -- the cross-pairing again.
    std::vector<std::string> partners = lz_;
    partners.insert(partners.end(), lx_.begin(), lx_.end());
    return class_of(sp_, partners, signature, "with_logical_signature");
}

// ===========================================================================
//  The change of coordinates
// ===========================================================================

Bdd StabilizerCode::to_code_coordinates(const Bdd &f) const
{
    Manager  &m      = sp_.manager();
    const int n_vars = sp_.n_vars();

    // The image of a set has characteristic function F(T^-1 y), so it is the
    // inverse that gets substituted: points push forward, functions pull back.
    // Each old variable becomes the XOR of the new ones its row selects, and
    // all 2n rules are applied at once because every right-hand side mentions
    // variables that are themselves being replaced.
    std::vector<std::pair<int, Bdd>> subs;
    subs.reserve(static_cast<std::size_t>(n_vars));
    for (int v = 0; v < n_vars; ++v) {
        Bdd g = m.constant(false);
        for (int rr = 0; rr < n_vars; ++rr)
            if (tinv_[static_cast<std::size_t>(v)][static_cast<std::size_t>(rr)])
                g ^= m.literal(rr, true);
        subs.emplace_back(v, std::move(g));
    }
    return f.compose(subs);
}

StabilizerCode::Row StabilizerCode::from_code_coordinates(const Row &y) const
{
    const int n_vars = sp_.n_vars();
    Row       e(static_cast<std::size_t>(n_vars), 0);
    for (int v = 0; v < n_vars; ++v)
        e[static_cast<std::size_t>(v)] = dot(tinv_[static_cast<std::size_t>(v)], y);
    return e;
}

// ===========================================================================
//  Inequivalent error pairs
// ===========================================================================
//
//  Two errors are indistinguishable to a decoder exactly when they share a
//  syndrome, and correcting one then leaves a logical error in the other case
//  exactly when their logical signatures differ. So the question is whether
//  any syndrome carries two signatures.
//
//  In the transformed coordinates the stabilizer component occupies its own
//  variables, so quotienting by S is just eliminating them, and what is left is
//  the (syndrome, signature) relation. Testing that relation for
//  single-valuedness would seem to need the 4^k signatures per syndrome; it
//  does not. Two distinct signatures must differ in *some* bit, so it is enough
//  to ask, one bit at a time, whether a syndrome admits both values of that
//  bit. That is 2k projections rather than 4^k.

namespace {

struct MultiTest {
    Bdd multi;      // the offending syndromes
    Bdd relation;   // the (syndrome, signature) set, for witness extraction
};

} // namespace

bool StabilizerCode::has_inequivalent_pair(const PauliSet &errors) const
{
    return find_inequivalent_pair(errors).has_value();
}

std::optional<StabilizerCode::Pair>
StabilizerCode::find_inequivalent_pair(const PauliSet &errors) const
{
    if (&errors.space().manager() != &sp_.manager())
        throw std::logic_error("find_inequivalent_pair: the set belongs to a different PauliSpace");
    if (k_ == 0) return std::nullopt;   // no logical operators, so no such pair
    if (errors.is_empty()) return std::nullopt;

    Manager &m = sp_.manager();

    std::vector<int> a_vars, lambda;
    for (int i = 0; i < r_; ++i) a_vars.push_back(a_var(i));
    for (int j = 0; j < k_; ++j) { lambda.push_back(c_var(j)); lambda.push_back(f_var(j)); }

    const Bdd transformed = to_code_coordinates(errors.bdd());
    const Bdd relation    = transformed.exists(a_vars);   // quotient by S

    // One bit at a time: the syndromes under which this signature bit takes
    // both values are syndromes carrying two distinct signatures.
    Bdd multi = m.constant(false);
    for (int v : lambda) {
        const Bdd zero = (relation & m.literal(v, false)).exists(lambda);
        const Bdd one  = (relation & m.literal(v, true)).exists(lambda);
        multi |= zero & one;
    }
    if (multi.is_false()) return std::nullopt;

    // --- witness extraction ------------------------------------------------
    // multi depends only on the syndrome variables, so any point of it names an
    // offending syndrome.
    auto bits_of = [&](const Bdd &f) {
        const PauliSet         one = sp_.wrap(f);
        const std::string      s   = *one.any_element();
        const std::vector<Pauli> p = parse_pauli_string(s);
        return to_row(p);
    };

    const Row offending = bits_of(multi);

    Bdd pinned_syndrome = m.constant(true);
    for (int i = 0; i < r_; ++i)
        pinned_syndrome &= m.literal(b_var(i), offending[static_cast<std::size_t>(b_var(i))] != 0);

    // Two signatures from that fiber: take one, then exclude it and take
    // another. The multi test guarantees a second exists.
    const Bdd fiber = relation & pinned_syndrome;
    const Row first_sig = bits_of(fiber);

    Bdd first_cube = m.constant(true);
    for (int v : lambda)
        first_cube &= m.literal(v, first_sig[static_cast<std::size_t>(v)] != 0);

    const Bdd rest = fiber.diff(first_cube);
    if (rest.is_false()) throw std::logic_error("find_inequivalent_pair: the multi test lied");
    const Row second_sig = bits_of(rest);

    // Back to actual operators: pin syndrome and signature in the transformed
    // set (before the quotient, so a stabilizer component is available), take
    // any point, and map it back through T^-1.
    auto witness = [&](const Row &sig) {
        Bdd cube = pinned_syndrome;
        for (int v : lambda) cube &= m.literal(v, sig[static_cast<std::size_t>(v)] != 0);
        const Row y = bits_of(transformed & cube);
        return row_to_string(from_code_coordinates(y));
    };

    Pair out;
    out.first            = witness(first_sig);
    out.second           = witness(second_sig);
    out.syndrome         = syndrome(out.first);
    out.signature_first  = logical_signature(out.first);
    out.signature_second = logical_signature(out.second);
    return out;
}

} // namespace spbdd
