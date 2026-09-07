// ===========================================================================
//  Implementation of spbdd::Bdd.
// ===========================================================================

#include "spbdd/bdd.hpp"
#include "spbdd/manager.hpp"

// cudd.h uses size_t and FILE but includes neither header itself, so these two
// lines are load-bearing and must come first.
#include <cstddef>
#include <cstdio>

#include <cudd.h>

#include <stdexcept>
#include <string>

namespace spbdd {
namespace {

// CUDD signals failure by returning null (typically out of memory, or a
// reordering that was interrupted). Turning that into an exception here keeps
// the null out of every call site.
DdNode *checked(DdNode *node, const char *what)
{
    if (!node) throw std::runtime_error(std::string("CUDD failed in ") + what);
    return node;
}

void require_valid(const Bdd &f, const char *what)
{
    if (!f.valid()) throw std::logic_error(std::string(what) + ": empty Bdd");
}

// Functions from two different Managers share no unique table, so combining
// them would produce a node that means nothing. Caught here rather than left
// to corrupt a run.
void require_same_manager(const Bdd &a, const Bdd &b, const char *what)
{
    require_valid(a, what);
    require_valid(b, what);
    if (a.manager() != b.manager())
        throw std::logic_error(std::string(what) + ": operands belong to different Managers");
}

void require_var_in_range(const Manager &mgr, int var, const char *what)
{
    if (var < 0 || var >= mgr.var_count())
        throw std::out_of_range(std::string(what) + ": variable " + std::to_string(var) +
                                " outside the declared space of " +
                                std::to_string(mgr.var_count()));
}

// The conjunction of the given variables, which is what CUDD's quantification
// calls take instead of a list.
Bdd make_cube(Manager *mgr, const std::vector<int> &vars, const char *what)
{
    for (int v : vars) require_var_in_range(*mgr, v, what);

    std::vector<int> indices(vars);   // Cudd_IndicesToCube wants a mutable int*
    DdNode          *cube = checked(Cudd_IndicesToCube(mgr->raw(), indices.data(),
                                                       static_cast<int>(indices.size())),
                                    "Cudd_IndicesToCube");
    Cudd_Ref(cube);
    return Bdd::adopt(mgr, cube);
}

} // namespace

// ===========================================================================
//  Lifetime
// ===========================================================================

Bdd Bdd::adopt(Manager *mgr, DdNode *node) noexcept
{
    return Bdd(mgr, node);
}

Bdd Bdd::borrow(Manager *mgr, DdNode *node)
{
    if (mgr && node) Cudd_Ref(node);
    return Bdd(mgr, node);
}

Bdd::Bdd(const Bdd &o) : mgr_(o.mgr_), node_(o.node_)
{
    if (mgr_) Cudd_Ref(node_);
}

Bdd::Bdd(Bdd &&o) noexcept : mgr_(o.mgr_), node_(o.node_)
{
    o.mgr_  = nullptr;
    o.node_ = nullptr;
}

Bdd &Bdd::operator=(const Bdd &o)
{
    if (this != &o) {
        // Reference first, release second: if both sides name the same node,
        // releasing first could drop its count to zero and free it.
        if (o.mgr_) Cudd_Ref(o.node_);
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
        o.node_ = nullptr;
    }
    return *this;
}

Bdd::~Bdd()
{
    release();
}

void Bdd::release() noexcept
{
    if (mgr_ && node_) Cudd_RecursiveDeref(mgr_->raw(), node_);
    mgr_  = nullptr;
    node_ = nullptr;
}

// ===========================================================================
//  Boolean algebra
// ===========================================================================

Bdd Bdd::operator!() const
{
    require_valid(*this, "operator!");
    DdNode *r = Cudd_Not(node_);
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd Bdd::operator&(const Bdd &o) const
{
    require_same_manager(*this, o, "operator&");
    DdNode *r = checked(Cudd_bddAnd(mgr_->raw(), node_, o.node_), "Cudd_bddAnd");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd Bdd::operator|(const Bdd &o) const
{
    require_same_manager(*this, o, "operator|");
    DdNode *r = checked(Cudd_bddOr(mgr_->raw(), node_, o.node_), "Cudd_bddOr");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd Bdd::operator^(const Bdd &o) const
{
    require_same_manager(*this, o, "operator^");
    DdNode *r = checked(Cudd_bddXor(mgr_->raw(), node_, o.node_), "Cudd_bddXor");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd Bdd::diff(const Bdd &o) const
{
    require_same_manager(*this, o, "diff");
    DdNode *r = checked(Cudd_bddAnd(mgr_->raw(), node_, Cudd_Not(o.node_)), "Cudd_bddAnd");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd Bdd::ite(const Bdd &then_case, const Bdd &else_case) const
{
    require_same_manager(*this, then_case, "ite");
    require_same_manager(*this, else_case, "ite");
    DdNode *r = checked(Cudd_bddIte(mgr_->raw(), node_, then_case.node_, else_case.node_),
                        "Cudd_bddIte");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd &Bdd::operator&=(const Bdd &o) { return *this = *this & o; }
Bdd &Bdd::operator|=(const Bdd &o) { return *this = *this | o; }
Bdd &Bdd::operator^=(const Bdd &o) { return *this = *this ^ o; }

// ===========================================================================
//  Quantification
// ===========================================================================

Bdd Bdd::exists(const std::vector<int> &vars) const
{
    require_valid(*this, "exists");
    if (vars.empty()) return *this;

    Bdd     cube = make_cube(mgr_, vars, "exists");
    DdNode *r    = checked(Cudd_bddExistAbstract(mgr_->raw(), node_, cube.node()),
                           "Cudd_bddExistAbstract");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd Bdd::forall(const std::vector<int> &vars) const
{
    require_valid(*this, "forall");
    if (vars.empty()) return *this;

    Bdd     cube = make_cube(mgr_, vars, "forall");
    DdNode *r    = checked(Cudd_bddUnivAbstract(mgr_->raw(), node_, cube.node()),
                           "Cudd_bddUnivAbstract");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

// ===========================================================================
//  Substitution
// ===========================================================================

Bdd Bdd::compose(const std::vector<std::pair<int, Bdd>> &subs) const
{
    require_valid(*this, "compose");
    if (subs.empty()) return *this;

    DdManager *dd = mgr_->raw();
    const int  n  = mgr_->var_count();

    // Cudd_bddVectorCompose wants an entry for *every* variable in the space,
    // not just the ones being replaced; anything left out has to be mapped to
    // itself explicitly.
    std::vector<DdNode *> vector(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        vector[static_cast<std::size_t>(i)] = checked(Cudd_bddIthVar(dd, i), "Cudd_bddIthVar");

    for (const auto &sub : subs) {
        require_var_in_range(*mgr_, sub.first, "compose");
        require_valid(sub.second, "compose");
        if (sub.second.manager() != mgr_)
            throw std::logic_error("compose: substituted function belongs to a different Manager");
        vector[static_cast<std::size_t>(sub.first)] = sub.second.node();
    }

    DdNode *r = checked(Cudd_bddVectorCompose(dd, node_, vector.data()), "Cudd_bddVectorCompose");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

Bdd Bdd::permute(const std::vector<std::pair<int, int>> &renaming) const
{
    require_valid(*this, "permute");
    if (renaming.empty()) return *this;

    const int n = mgr_->var_count();

    // Same rule as compose: the array covers the whole space, so start from the
    // identity and overwrite only what the caller named.
    std::vector<int> permutation(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) permutation[static_cast<std::size_t>(i)] = i;

    for (const auto &rename : renaming) {
        require_var_in_range(*mgr_, rename.first, "permute");
        require_var_in_range(*mgr_, rename.second, "permute");
        permutation[static_cast<std::size_t>(rename.first)] = rename.second;
    }

    DdNode *r = checked(Cudd_bddPermute(mgr_->raw(), node_, permutation.data()),
                        "Cudd_bddPermute");
    Cudd_Ref(r);
    return Bdd(mgr_, r);
}

// ===========================================================================
//  Structure
// ===========================================================================

bool Bdd::is_constant() const
{
    require_valid(*this, "is_constant");
    return Cudd_IsConstant(node_) != 0;
}

bool Bdd::is_true() const
{
    require_valid(*this, "is_true");
    return node_ == Cudd_ReadOne(mgr_->raw());
}

bool Bdd::is_false() const
{
    require_valid(*this, "is_false");
    return node_ == Cudd_ReadLogicZero(mgr_->raw());
}

int Bdd::top_var() const
{
    require_valid(*this, "top_var");
    if (Cudd_IsConstant(node_)) throw std::logic_error("top_var: node is a constant");
    return static_cast<int>(Cudd_NodeReadIndex(node_));
}

Bdd Bdd::low() const
{
    require_valid(*this, "low");
    if (Cudd_IsConstant(node_)) throw std::logic_error("low: node is a constant");

    // Cudd_E strips the mark off the parent edge but does not push it down, so
    // a complemented parent has to hand a complemented child back.
    DdNode *e = Cudd_E(node_);
    if (Cudd_IsComplement(node_)) e = Cudd_Not(e);
    Cudd_Ref(e);
    return Bdd(mgr_, e);
}

Bdd Bdd::high() const
{
    require_valid(*this, "high");
    if (Cudd_IsConstant(node_)) throw std::logic_error("high: node is a constant");

    DdNode *t = Cudd_T(node_);
    if (Cudd_IsComplement(node_)) t = Cudd_Not(t);
    Cudd_Ref(t);
    return Bdd(mgr_, t);
}

// ===========================================================================
//  Measurement
// ===========================================================================

double Bdd::sat_count(int n_vars) const
{
    require_valid(*this, "sat_count");
    if (n_vars < 0) throw std::invalid_argument("sat_count: negative variable count");
    return Cudd_CountMinterm(mgr_->raw(), node_, n_vars);
}

std::size_t Bdd::node_count() const
{
    require_valid(*this, "node_count");
    return static_cast<std::size_t>(Cudd_DagSize(node_));
}

} // namespace spbdd
