// ===========================================================================
//  PauliSet -- set algebra, predicates, metrics and enumeration.
//
//  Everything except the enumeration block is a forwarding call to Bdd: the
//  work was already done by the builders in paulispace.cpp, and a set
//  operation on characteristic functions *is* the boolean operation.
// ===========================================================================

#include "spbdd/pauliset.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace spbdd {
namespace {

void require_qubit(const PauliSpace &sp, int q, const char *who)
{
    if (q < 0 || q >= sp.n_qubits())
        throw std::out_of_range(std::string(who) + ": qubit " + std::to_string(q) +
                                " outside a space of " + std::to_string(sp.n_qubits()));
}

void require_same_space(const PauliSet &a, const PauliSet &b, const char *who)
{
    if (&a.space().manager() != &b.space().manager())
        throw std::logic_error(std::string(who) + ": operands come from different PauliSpaces");
}

} // namespace

// ===========================================================================
//  Set algebra
// ===========================================================================

PauliSet PauliSet::operator|(const PauliSet &o) const
{
    require_same_space(*this, o, "operator|");
    return sp_.wrap(f_ | o.f_);
}

PauliSet PauliSet::operator&(const PauliSet &o) const
{
    require_same_space(*this, o, "operator&");
    return sp_.wrap(f_ & o.f_);
}

PauliSet PauliSet::operator-(const PauliSet &o) const
{
    require_same_space(*this, o, "operator-");
    return sp_.wrap(f_.diff(o.f_));
}

PauliSet PauliSet::operator^(const PauliSet &o) const
{
    require_same_space(*this, o, "operator^");
    return sp_.wrap(f_ ^ o.f_);
}

PauliSet PauliSet::operator~() const { return sp_.wrap(!f_); }

PauliSet &PauliSet::operator|=(const PauliSet &o) { return *this = *this | o; }
PauliSet &PauliSet::operator&=(const PauliSet &o) { return *this = *this & o; }
PauliSet &PauliSet::operator-=(const PauliSet &o) { return *this = *this - o; }
PauliSet &PauliSet::operator^=(const PauliSet &o) { return *this = *this ^ o; }

bool PauliSet::operator==(const PauliSet &o) const
{
    return &sp_.manager() == &o.sp_.manager() && f_ == o.f_;
}

// ===========================================================================
//  Predicates
// ===========================================================================

bool PauliSet::is_empty() const { return f_.is_false(); }

// The universe is the constant true because there is no unused encoding.
bool PauliSet::is_universe() const { return f_.is_true(); }

bool PauliSet::contains(const std::vector<Pauli> &p) const
{
    return !(f_ & sp_.from(p).bdd()).is_false();
}

bool PauliSet::contains(const std::string &p) const
{
    return contains(parse_pauli_string(p, sp_.n_qubits(), "contains"));
}

bool PauliSet::subset_of(const PauliSet &o) const
{
    require_same_space(*this, o, "subset_of");
    return f_.diff(o.f_).is_false();
}

bool PauliSet::superset_of(const PauliSet &o) const { return o.subset_of(*this); }

bool PauliSet::disjoint_from(const PauliSet &o) const
{
    require_same_space(*this, o, "disjoint_from");
    return (f_ & o.f_).is_false();
}

// ===========================================================================
//  Metrics
// ===========================================================================

double      PauliSet::size() const { return f_.sat_count(sp_.n_vars()); }
std::size_t PauliSet::node_count() const { return f_.node_count(); }

int PauliSet::min_weight() const
{
    if (is_empty()) return -1;

    // "some element has weight <= w" is monotone in w and true at w = n, so
    // the smallest weight present is found by bisection. Each probe builds one
    // weight counter; nothing is ever enumerated.
    int lo = 0;
    int hi = sp_.n_qubits();
    while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (!(f_ & sp_.weight_at_most(mid).bdd()).is_false()) hi = mid;
        else                                                  lo = mid + 1;
    }
    return lo;
}

// ===========================================================================
//  Derived sets
// ===========================================================================

PauliSet PauliSet::forget(const std::vector<int> &qubits) const
{
    std::vector<int> vars;
    vars.reserve(2 * qubits.size());
    for (int q : qubits) {
        require_qubit(sp_, q, "forget");
        vars.push_back(PauliSpace::xvar(q));
        vars.push_back(PauliSpace::zvar(q));
    }
    return sp_.wrap(f_.exists(vars));
}

// A fault at a set of locations can leave behind any Pauli supported on
// them, which is exactly what forget() does. Named separately because this
// is what a circuit-level caller is looking for.
PauliSet PauliSet::fault_inject(const std::vector<int> &qubits) const
{
    return forget(qubits);
}

PauliSet PauliSet::reset(const std::vector<int> &qubits) const
{
    // Quantifying guarantees that the "identity on these qubits" member of
    // every coset is present; intersecting then selects it. Each coset has
    // exactly one such member, so nothing is lost or duplicated.
    return sp_.wrap(forget(qubits).bdd() & sp_.identity_on(qubits).bdd());
}

PauliSet PauliSet::with_weight_at_most(int w) const
{
    return sp_.wrap(f_ & sp_.weight_at_most(w).bdd());
}

PauliSet PauliSet::with_weight_exactly(int w) const
{
    return sp_.wrap(f_ & sp_.weight_exactly(w).bdd());
}

// ===========================================================================
//  Clifford gates
// ===========================================================================
//
//  A Clifford conjugation is determined by where it sends the generators
//  X_q and Z_q, and it extends linearly from there. Writing each image in the
//  (x, z) coordinates gives a linear map on GF(2)^(2n), which is exactly what
//  Bdd::compose and Bdd::permute apply -- in one simultaneous step, which for
//  a map whose right-hand sides mention each other is a correctness
//  requirement rather than an optimisation.
//
//  Substituting the map itself computes f . M, while the image of a set under
//  M is f . M^-1. The two agree here because every gate below squares to a
//  Pauli, which acts trivially in the phase-free representation, so M is an
//  involution. A gate without that property would have to substitute the
//  inverse map.

namespace {

void require_two_qubits(const PauliSpace &sp, int a, int b, const char *who)
{
    require_qubit(sp, a, who);
    require_qubit(sp, b, who);
    if (a == b)
        throw std::invalid_argument(std::string(who) + ": the two qubits must differ");
}

} // namespace

// Conjugation by a Pauli changes only a sign, which this representation does
// not carry. The qubit is still checked, so a circuit with a bad index fails
// here rather than silently doing nothing.
PauliSet PauliSet::x(int qubit) const
{
    require_qubit(sp_, qubit, "x");
    return *this;
}

PauliSet PauliSet::y(int qubit) const
{
    require_qubit(sp_, qubit, "y");
    return *this;
}

PauliSet PauliSet::z(int qubit) const
{
    require_qubit(sp_, qubit, "z");
    return *this;
}

// X -> Z, Z -> X: the two coordinates of the qubit trade places. A renaming,
// so permute rather than compose.
PauliSet PauliSet::h(int qubit) const
{
    require_qubit(sp_, qubit, "h");
    const int xv = PauliSpace::xvar(qubit);
    const int zv = PauliSpace::zvar(qubit);
    return sp_.wrap(f_.permute({{xv, zv}, {zv, xv}}));
}

// X -> Y = XZ, Z -> Z. In coordinates z_q gains x_q.
PauliSet PauliSet::s(int qubit) const
{
    require_qubit(sp_, qubit, "s");
    Manager  &m  = sp_.manager();
    const int xv = PauliSpace::xvar(qubit);
    const int zv = PauliSpace::zvar(qubit);
    return sp_.wrap(f_.compose({{zv, m.literal(xv) ^ m.literal(zv)}}));
}

// S* sends X -> Y, Z -> Z as well: S and S* differ only by the sign of Y,
// which this phase-free representation does not carry. A separate entry point
// exists so a circuit can be replayed without special-casing adjoints.
PauliSet PauliSet::sdg(int qubit) const
{
    require_qubit(sp_, qubit, "sdg");
    Manager  &m  = sp_.manager();
    const int xv = PauliSpace::xvar(qubit);
    const int zv = PauliSpace::zvar(qubit);
    return sp_.wrap(f_.compose({{zv, m.literal(xv) ^ m.literal(zv)}}));
}

// X -> X, Z -> Y = XZ. The mirror image of s().
PauliSet PauliSet::sx(int qubit) const
{
    require_qubit(sp_, qubit, "sx");
    Manager  &m  = sp_.manager();
    const int xv = PauliSpace::xvar(qubit);
    const int zv = PauliSpace::zvar(qubit);
    return sp_.wrap(f_.compose({{xv, m.literal(xv) ^ m.literal(zv)}}));
}

// sqrt(X)* sends X -> X, Z -> Y as well, for the same reason sdg() agrees
// with s(): the sign that distinguishes the adjoint is not tracked here.
PauliSet PauliSet::sxdg(int qubit) const
{
    require_qubit(sp_, qubit, "sxdg");
    Manager  &m  = sp_.manager();
    const int xv = PauliSpace::xvar(qubit);
    const int zv = PauliSpace::zvar(qubit);
    return sp_.wrap(f_.compose({{xv, m.literal(xv) ^ m.literal(zv)}}));
}

// X_c -> X_c X_t,  Z_c -> Z_c,  X_t -> X_t,  Z_t -> Z_c Z_t.
// So x flows from control to target and z flows back the other way -- the two
// assignments read each other, which is why they must happen together.
PauliSet PauliSet::cx(int control, int target) const
{
    require_two_qubits(sp_, control, target, "cx");
    Manager  &m  = sp_.manager();
    const int xc = PauliSpace::xvar(control);
    const int zc = PauliSpace::zvar(control);
    const int xt = PauliSpace::xvar(target);
    const int zt = PauliSpace::zvar(target);

    return sp_.wrap(f_.compose({{xt, m.literal(xc) ^ m.literal(xt)},
                                {zc, m.literal(zc) ^ m.literal(zt)}}));
}

// X_c -> X_c Y_t,  Z_c -> Z_c,  X_t -> Z_c X_t,  Z_t -> Z_c Z_t.
PauliSet PauliSet::cy(int control, int target) const
{
    require_two_qubits(sp_, control, target, "cy");
    Manager  &m  = sp_.manager();
    const int xc = PauliSpace::xvar(control);
    const int zc = PauliSpace::zvar(control);
    const int xt = PauliSpace::xvar(target);
    const int zt = PauliSpace::zvar(target);

    return sp_.wrap(f_.compose({{zc, m.literal(zc) ^ m.literal(xt) ^ m.literal(zt)},
                                {xt, m.literal(xc) ^ m.literal(xt)},
                                {zt, m.literal(xc) ^ m.literal(zt)}}));
}

// X_a -> X_a Z_b,  X_b -> Z_a X_b,  both Z generators fixed. Symmetric in its
// two arguments, unlike cx and cy.
PauliSet PauliSet::cz(int a, int b) const
{
    require_two_qubits(sp_, a, b, "cz");
    Manager  &m  = sp_.manager();
    const int xa = PauliSpace::xvar(a);
    const int za = PauliSpace::zvar(a);
    const int xb = PauliSpace::xvar(b);
    const int zb = PauliSpace::zvar(b);

    return sp_.wrap(f_.compose({{za, m.literal(za) ^ m.literal(xb)},
                                {zb, m.literal(zb) ^ m.literal(xa)}}));
}

PauliSet PauliSet::swap(int a, int b) const
{
    require_two_qubits(sp_, a, b, "swap");
    const int xa = PauliSpace::xvar(a);
    const int za = PauliSpace::zvar(a);
    const int xb = PauliSpace::xvar(b);
    const int zb = PauliSpace::zvar(b);

    return sp_.wrap(f_.permute({{xa, xb}, {za, zb}, {xb, xa}, {zb, za}}));
}



// ===========================================================================
//  Measurement
// ===========================================================================

namespace {

// The observable as an n-character Pauli string, for the one-qubit cases.
std::string single_qubit_observable(const PauliSpace &sp, int qubit, char letter, const char *who)
{
    require_qubit(sp, qubit, who);
    std::string s(static_cast<std::size_t>(sp.n_qubits()), 'I');
    s[static_cast<std::size_t>(qubit)] = letter;
    return s;
}

} // namespace

MeasurementSplit PauliSet::measure(const std::string &observable) const
{
    // "Anticommutes with the observable" is one parity check, so each side is a
    // single intersection and neither one enumerates anything. Together the two
    // sides partition the set: an empty side is an outcome that cannot occur.
    return MeasurementSplit{*this & sp_.commuting_with(observable),
                            *this & sp_.anticommuting_with(observable)};
}

// X and Y anticommute with Z while I and Z do not, so a Z-basis measurement
// reads the x coordinate of the qubit.
MeasurementSplit PauliSet::measure_z(int qubit) const
{
    return measure(single_qubit_observable(sp_, qubit, 'Z', "measure_z"));
}

MeasurementSplit PauliSet::measure_x(int qubit) const
{
    return measure(single_qubit_observable(sp_, qubit, 'X', "measure_x"));
}

MeasurementSplit PauliSet::measure_y(int qubit) const
{
    return measure(single_qubit_observable(sp_, qubit, 'Y', "measure_y"));
}

// ===========================================================================
//  Enumeration
// ===========================================================================
//
//  A path through the diagram skips every variable the function does not
//  depend on, so one path stands for 2^k assignments. Expanding those skipped
//  variables here is what makes each callback a single concrete operator
//  rather than a cube with don't-cares in it.
//
//  The walk is by level, not by variable number, because that is the order the
//  diagram is built in; which is also why the order of the callbacks is not
//  promised and to_strings() sorts.

namespace {

class Enumerator {
public:
    Enumerator(Manager &mgr, int n_qubits,
               const std::function<bool(const std::vector<Pauli> &)> &fn)
        : mgr_(mgr), n_qubits_(n_qubits), n_vars_(2 * n_qubits),
          bits_(static_cast<std::size_t>(2 * n_qubits), 0),
          element_(static_cast<std::size_t>(n_qubits), Pauli::I), fn_(fn)
    {}

    // Returns false once the caller has asked to stop.
    bool walk(const Bdd &f, int level)
    {
        if (f.is_false()) return true;

        if (level == n_vars_) {
            for (int q = 0; q < n_qubits_; ++q)
                element_[static_cast<std::size_t>(q)] =
                    make_pauli(bits_[static_cast<std::size_t>(2 * q)] != 0,
                               bits_[static_cast<std::size_t>(2 * q + 1)] != 0);
            return fn_(element_);
        }

        const int  var     = mgr_.level_to_var(level);
        const bool asked   = !f.is_constant() && mgr_.var_to_level(f.top_var()) == level;
        const auto recurse = [&](char value, const Bdd &child) {
            bits_[static_cast<std::size_t>(var)] = value;
            return walk(child, level + 1);
        };

        if (!asked)   // a don't-care: both values are in the set
            return recurse(0, f) && recurse(1, f);

        return recurse(0, f.low()) && recurse(1, f.high());
    }

private:
    Manager           &mgr_;
    int                n_qubits_;
    int                n_vars_;
    std::vector<char>  bits_;      // indexed by variable number
    std::vector<Pauli> element_;

    const std::function<bool(const std::vector<Pauli> &)> &fn_;
};

} // namespace

void PauliSet::for_each(const std::function<bool(const std::vector<Pauli> &)> &fn) const
{
    Enumerator e(sp_.manager(), sp_.n_qubits(), fn);
    e.walk(f_, 0);
}

// The same walk, handing each element over as a string. Most callers want
// this one; the vector<Pauli> form avoids the formatting when they do not.
void PauliSet::for_each(const std::function<bool(const std::string &)> &fn) const
{
    for_each([&](const std::vector<Pauli> &p) { return fn(pauli_string_to_text(p)); });
}

std::vector<std::string> PauliSet::to_strings(std::size_t max_elements) const
{
    std::vector<std::string> out;
    bool                     overflowed = false;

    for_each([&](const std::vector<Pauli> &p) {
        if (out.size() >= max_elements) {
            overflowed = true;
            return false;
        }
        out.push_back(pauli_string_to_text(p));
        return true;
    });

    if (overflowed)
        throw std::length_error("to_strings: more than " + std::to_string(max_elements) +
                                " elements");

    std::sort(out.begin(), out.end());
    return out;
}

std::optional<std::string> PauliSet::any_element() const
{
    if (is_empty()) return std::nullopt;

    Manager  &mgr    = sp_.manager();
    const int n_vars = sp_.n_vars();

    std::vector<char> bits(static_cast<std::size_t>(n_vars), 0);
    Bdd               cur = f_;

    // One pass down the diagram. Where the function does not ask about a
    // variable, either value will do; where it does, take whichever branch is
    // still satisfiable.
    for (int level = 0; level < n_vars; ++level) {
        const int var = mgr.level_to_var(level);
        if (cur.is_constant() || mgr.var_to_level(cur.top_var()) != level) continue;

        Bdd low = cur.low();
        if (!low.is_false()) {
            bits[static_cast<std::size_t>(var)] = 0;
            cur                                 = std::move(low);
        } else {
            bits[static_cast<std::size_t>(var)] = 1;
            cur                                 = cur.high();
        }
    }

    std::vector<Pauli> p(static_cast<std::size_t>(sp_.n_qubits()), Pauli::I);
    for (int q = 0; q < sp_.n_qubits(); ++q)
        p[static_cast<std::size_t>(q)] = make_pauli(bits[static_cast<std::size_t>(2 * q)] != 0,
                                                    bits[static_cast<std::size_t>(2 * q + 1)] != 0);
    return pauli_string_to_text(p);
}

} // namespace spbdd
