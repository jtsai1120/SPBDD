// ===========================================================================
//  Implementation of spbdd::Bdd against BuDDy.
//
//  BuDDy hands back an unreferenced node and expects the caller to claim it
//  before anything else can trigger a garbage collection, so every result here
//  is bdd_addref'd on the spot and every intermediate is released.
//
//  BuDDy's error handler prints and terminates the process on an unrecoverable
//  condition (node table exhausted, illegal variable). That behaviour is left
//  alone: throwing a C++ exception back through BuDDy's C frames would leave
//  its tables inconsistent and would be worse than stopping.
// ===========================================================================

#include "spbdd/bdd.hpp"
#include "spbdd/manager.hpp"

#include <bdd.h>

// BuDDy's header, compiled as C++, redirects part of its C API to overloads
// returning objects of its own bdd class -- including the terminals and the
// variable constructors this layer uses. The C entry points are still in the
// library, so take them back: everything below wants plain node indices.
#undef bddtrue
#undef bddfalse
#undef bdd_init
#undef bdd_ithvar
#undef bdd_nithvar
#undef bdd_makeset
#undef bdd_ibuildcube
#undef bdd_anodecount
extern "C" const BDD bddtrue;
extern "C" const BDD bddfalse;

#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace spbdd {
namespace {

NodeId claim(NodeId node)
{
    bdd_addref(node);
    return node;
}

void require_valid(const Bdd &f, const char *what)
{
    if (!f.valid()) throw std::logic_error(std::string(what) + ": empty Bdd");
}

// With BuDDy there is only ever one manager, so this reduces to "neither side
// is empty" -- but the check is kept so the contract reads the same as the
// CUDD backend's, and so an empty Bdd is caught rather than treated as false.
void require_same_manager(const Bdd &a, const Bdd &b, const char *what)
{
    require_valid(a, what);
    require_valid(b, what);
    if (a.manager() != b.manager())
        throw std::logic_error(std::string(what) + ": operands belong to different Managers");
}

void require_var_in_range(int var, const char *what)
{
    if (var < 0 || var >= bdd_varnum())
        throw std::out_of_range(std::string(what) + ": variable " + std::to_string(var) +
                                " outside the declared space of " +
                                std::to_string(bdd_varnum()));
}

// The conjunction of the given variables, which is what BuDDy's quantification
// calls take instead of a list. Built explicitly rather than with bdd_makeset,
// which wants its input in level order.
Bdd make_cube(Manager *mgr, const std::vector<int> &vars, const char *what)
{
    NodeId cube = claim(bddtrue);
    for (int v : vars) {
        require_var_in_range(v, what);
        const NodeId next = claim(bdd_and(cube, bdd_ithvar(v)));
        bdd_delref(cube);
        cube = next;
    }
    return Bdd::adopt(mgr, cube);
}

} // namespace

// ===========================================================================
//  Lifetime
// ===========================================================================

Bdd Bdd::adopt(Manager *mgr, NodeId node) noexcept { return Bdd(mgr, node); }

Bdd Bdd::borrow(Manager *mgr, NodeId node)
{
    if (mgr) bdd_addref(node);
    return Bdd(mgr, node);
}

Bdd::Bdd(const Bdd &o) : mgr_(o.mgr_), node_(o.node_)
{
    if (mgr_) bdd_addref(node_);
}

Bdd::Bdd(Bdd &&o) noexcept : mgr_(o.mgr_), node_(o.node_)
{
    o.mgr_  = nullptr;
    o.node_ = 0;
}

Bdd &Bdd::operator=(const Bdd &o)
{
    if (this != &o) {
        // Reference first, release second: if both sides name the same node,
        // releasing first could drop its count to zero and free it.
        if (o.mgr_) bdd_addref(o.node_);
        release();
        mgr_  = o.mgr_;
        node_ = o.node_;
    }
    return *this;
}

Bdd &Bdd::operator=(Bdd &&o) noexcept
{
    if (this != &o) {
        release();
        mgr_    = o.mgr_;
        node_   = o.node_;
        o.mgr_  = nullptr;
        o.node_ = 0;
    }
    return *this;
}

Bdd::~Bdd() { release(); }

void Bdd::release() noexcept
{
    // Only touch BuDDy while a Manager is alive: a Bdd that outlives its
    // Manager would otherwise call into a torn-down node table.
    if (mgr_ && Manager::alive()) bdd_delref(node_);
    mgr_  = nullptr;
    node_ = 0;
}

// ===========================================================================
//  Boolean algebra
// ===========================================================================

Bdd Bdd::operator!() const
{
    require_valid(*this, "operator!");
    return Bdd(mgr_, claim(bdd_not(node_)));
}

Bdd Bdd::operator&(const Bdd &o) const
{
    require_same_manager(*this, o, "operator&");
    return Bdd(mgr_, claim(bdd_and(node_, o.node_)));
}

Bdd Bdd::operator|(const Bdd &o) const
{
    require_same_manager(*this, o, "operator|");
    return Bdd(mgr_, claim(bdd_or(node_, o.node_)));
}

Bdd Bdd::operator^(const Bdd &o) const
{
    require_same_manager(*this, o, "operator^");
    return Bdd(mgr_, claim(bdd_xor(node_, o.node_)));
}

Bdd Bdd::diff(const Bdd &o) const
{
    require_same_manager(*this, o, "diff");
    return Bdd(mgr_, claim(bdd_apply(node_, o.node_, bddop_diff)));
}

Bdd Bdd::ite(const Bdd &then_case, const Bdd &else_case) const
{
    require_same_manager(*this, then_case, "ite");
    require_same_manager(*this, else_case, "ite");
    return Bdd(mgr_, claim(bdd_ite(node_, then_case.node_, else_case.node_)));
}

Bdd &Bdd::operator&=(const Bdd &o) { return *this = *this & o; }
Bdd &Bdd::operator|=(const Bdd &o) { return *this = *this | o; }
Bdd &Bdd::operator^=(const Bdd &o) { return *this = *this ^ o; }

// ===========================================================================
//  Sumset
// ===========================================================================
//
//  XOR acts on every coordinate independently, so the top variable v can be
//  peeled off: a pair (a, b) lands on the v=1 side of the result exactly when
//  a and b disagree on v, and on the v=0 side when they agree. Hence
//
//      F * G = ITE(v,  F0*G1 | F1*G0,  F0*G0 | F1*G1)
//
//  which is apply's skeleton with all four cofactor pairs instead of the two
//  diagonal ones. The only BDD operations are or and ite; the xor is entirely
//  in which cofactor feeds which side. When only one operand depends on v the
//  two sides coincide and the ite is skipped.
//
//  Results are memoised on the unordered node pair -- canonicity makes that
//  sound, commutativity makes it unordered -- so there are at most |F|*|G|
//  subproblems. Levels have to hold still for the whole recursion, so dynamic
//  reordering is switched off around it.

namespace {

struct NodePairHash {
    std::size_t operator()(const std::pair<int, int> &p) const noexcept
    {
        return std::hash<long long>()((static_cast<long long>(p.first) << 32) |
                                      static_cast<unsigned int>(p.second));
    }
};

class Sumset {
public:
    explicit Sumset(Manager *mgr) : mgr_(mgr) {}

    Bdd run(const Bdd &f, const Bdd &g)
    {
        if (f.is_false() || g.is_false()) return mgr_->constant(false);
        // The other operand is nonempty here, and the whole space times
        // anything nonempty is the whole space.
        if (f.is_true() || g.is_true()) return mgr_->constant(true);

        const std::pair<int, int> key = f.node() < g.node()
                                            ? std::make_pair(f.node(), g.node())
                                            : std::make_pair(g.node(), f.node());
        const auto hit = memo_.find(key);
        if (hit != memo_.end()) return hit->second;

        const int lf = bdd_var2level(f.top_var());
        const int lg = bdd_var2level(g.top_var());

        Bdd result;
        if (lf < lg) {
            result = run(f.low(), g) | run(f.high(), g);
        } else if (lg < lf) {
            result = run(f, g.low()) | run(f, g.high());
        } else {
            const Bdd f0 = f.low(), f1 = f.high();
            const Bdd g0 = g.low(), g1 = g.high();
            const Bdd same  = run(f0, g0) | run(f1, g1);
            const Bdd cross = run(f0, g1) | run(f1, g0);
            result = mgr_->literal(f.top_var()).ite(cross, same);
        }
        memo_.emplace(key, result);
        return result;
    }

private:
    Manager                                                   *mgr_;
    std::unordered_map<std::pair<int, int>, Bdd, NodePairHash> memo_;
};

} // namespace

Bdd Bdd::sumset(const Bdd &o) const
{
    require_same_manager(*this, o, "sumset");

    bdd_disable_reorder();
    Bdd result;
    try {
        result = Sumset(mgr_).run(*this, o);
    } catch (...) {
        bdd_enable_reorder();
        throw;
    }
    bdd_enable_reorder();
    return result;
}

// ===========================================================================
//  Quantification
// ===========================================================================

Bdd Bdd::exists(const std::vector<int> &vars) const
{
    require_valid(*this, "exists");
    if (vars.empty()) return *this;

    const Bdd cube = make_cube(mgr_, vars, "exists");
    return Bdd(mgr_, claim(bdd_exist(node_, cube.node())));
}

Bdd Bdd::forall(const std::vector<int> &vars) const
{
    require_valid(*this, "forall");
    if (vars.empty()) return *this;

    const Bdd cube = make_cube(mgr_, vars, "forall");
    return Bdd(mgr_, claim(bdd_forall(node_, cube.node())));
}

// ===========================================================================
//  Substitution
// ===========================================================================

Bdd Bdd::compose(const std::vector<std::pair<int, Bdd>> &subs) const
{
    require_valid(*this, "compose");
    if (subs.empty()) return *this;

    // A fresh pair starts as the identity on every variable, so unlike CUDD's
    // vector compose only the entries being replaced have to be filled in.
    bddPair *pair = bdd_newpair();
    if (!pair) throw std::runtime_error("bdd_newpair failed");

    try {
        for (const auto &sub : subs) {
            require_var_in_range(sub.first, "compose");
            require_valid(sub.second, "compose");
            if (sub.second.manager() != mgr_)
                throw std::logic_error(
                    "compose: substituted function belongs to a different Manager");
            bdd_setbddpair(pair, sub.first, sub.second.node_);
        }
        const NodeId r = claim(bdd_veccompose(node_, pair));
        bdd_freepair(pair);
        return Bdd(mgr_, r);
    } catch (...) {
        bdd_freepair(pair);
        throw;
    }
}

Bdd Bdd::permute(const std::vector<std::pair<int, int>> &renaming) const
{
    require_valid(*this, "permute");
    if (renaming.empty()) return *this;

    bddPair *pair = bdd_newpair();
    if (!pair) throw std::runtime_error("bdd_newpair failed");

    try {
        // Simultaneous, so a transposition comes out right rather than
        // collapsing both variables onto one.
        for (const auto &rename : renaming) {
            require_var_in_range(rename.first, "permute");
            require_var_in_range(rename.second, "permute");
            bdd_setpair(pair, rename.first, rename.second);
        }
        const NodeId r = claim(bdd_replace(node_, pair));
        bdd_freepair(pair);
        return Bdd(mgr_, r);
    } catch (...) {
        bdd_freepair(pair);
        throw;
    }
}

// ===========================================================================
//  Structure
// ===========================================================================

bool Bdd::is_constant() const
{
    require_valid(*this, "is_constant");
    return node_ == bddtrue || node_ == bddfalse;
}

bool Bdd::is_true() const
{
    require_valid(*this, "is_true");
    return node_ == bddtrue;
}

bool Bdd::is_false() const
{
    require_valid(*this, "is_false");
    return node_ == bddfalse;
}

int Bdd::top_var() const
{
    require_valid(*this, "top_var");
    if (is_constant()) throw std::logic_error("top_var: node is a constant");
    return bdd_var(node_);
}

// BuDDy has no complement edges, so the children are already the plain Shannon
// decomposition with no mark to push down.
Bdd Bdd::low() const
{
    require_valid(*this, "low");
    if (is_constant()) throw std::logic_error("low: node is a constant");
    return Bdd(mgr_, claim(bdd_low(node_)));
}

Bdd Bdd::high() const
{
    require_valid(*this, "high");
    if (is_constant()) throw std::logic_error("high: node is a constant");
    return Bdd(mgr_, claim(bdd_high(node_)));
}

// ===========================================================================
//  Measurement
// ===========================================================================

double Bdd::sat_count(int n_vars) const
{
    require_valid(*this, "sat_count");
    if (n_vars < 0) throw std::invalid_argument("sat_count: negative variable count");

    // BuDDy counts over every declared variable rather than over a space the
    // caller names, so the answer is rescaled. That is exact whenever the
    // function does not depend on the variables outside n_vars, which is the
    // only way this library uses it.
    const double count    = bdd_satcount(node_);
    const int    declared = bdd_varnum();
    if (n_vars == declared) return count;
    return count * std::pow(2.0, n_vars - declared);
}

std::size_t Bdd::node_count() const
{
    require_valid(*this, "node_count");
    // BuDDy counts internal nodes only; the terminal is counted here so that a
    // constant is 1 and "x0 and x1" is 3, as under CUDD.
    return static_cast<std::size_t>(bdd_nodecount(node_)) + 1u;
}

} // namespace spbdd
