#pragma once

// ===========================================================================
//  spbdd::Bdd -- an owning handle on one BDD node.
//
//  CUDD's own header is deliberately absent here. DdNode is forward-declared,
//  so a program that uses spbdd needs -Iinclude and nothing else on its
//  include path. The price is that every function which has to call into CUDD
//  is defined in src/bdd.cpp instead of inline; next to the cost of a BDD
//  operation that is not measurable.
//
//  Two properties everything above this layer relies on:
//
//    1. Reference counting. A Bdd owns exactly one reference to its node and
//       hands it back when it dies. No code above this file writes a
//       Cudd_Ref by hand.
//
//    2. Canonicity. Two Bdds are equal iff they denote the same boolean
//       function, so equality is a pointer comparison, never a traversal.
//
//  A Bdd points at its Manager without keeping it alive; the Manager must
//  outlive every Bdd made from it. That is not a burden in practice because a
//  Bdd is always reached through a PauliSet, which does own its space.
//
//  Nothing here can be constructed on its own -- see manager.hpp for the
//  factory, or include spbdd.hpp and forget the distinction.
// ===========================================================================

#include <cstddef>
#include <utility>
#include <vector>

struct DdNode;      // CUDD's node type, forward-declared on purpose (see above)

namespace spbdd {

class Manager;

class Bdd {
public:
    Bdd() noexcept = default;

    // Takes over a reference the caller already owns -- that is, exactly what
    // a CUDD call returns once it has been Cudd_Ref'd.
    static Bdd adopt(Manager *mgr, DdNode *node) noexcept;

    // Adds a reference of its own.
    static Bdd borrow(Manager *mgr, DdNode *node);

    Bdd(const Bdd &);
    Bdd(Bdd &&) noexcept;
    Bdd &operator=(const Bdd &);
    Bdd &operator=(Bdd &&) noexcept;
    ~Bdd();

    bool     valid() const noexcept { return mgr_ != nullptr; }
    Manager *manager() const noexcept { return mgr_; }
    DdNode  *node() const noexcept { return node_; }

    // Canonicity (see the header comment) makes this a pointer comparison.
    bool operator==(const Bdd &o) const noexcept { return mgr_ == o.mgr_ && node_ == o.node_; }
    bool operator!=(const Bdd &o) const noexcept { return !(*this == o); }

    // --- boolean algebra ---------------------------------------------------
    Bdd operator!() const;
    Bdd operator&(const Bdd &o) const;
    Bdd operator|(const Bdd &o) const;
    Bdd operator^(const Bdd &o) const;
    Bdd diff(const Bdd &o) const;                        // *this AND NOT o
    Bdd ite(const Bdd &then_case, const Bdd &else_case) const;

    Bdd &operator&=(const Bdd &o);
    Bdd &operator|=(const Bdd &o);
    Bdd &operator^=(const Bdd &o);

    // --- quantification ----------------------------------------------------
    // `vars` is a plain list of variable numbers; the cube CUDD wants is built
    // and released internally.
    Bdd exists(const std::vector<int> &vars) const;
    Bdd forall(const std::vector<int> &vars) const;

    // --- substitution ------------------------------------------------------
    // Both are *simultaneous*: every right-hand side is read in the original
    // function, never in a partially substituted one. That is what makes a map
    // like {a := a xor b, b := a} come out right, and it is a correctness
    // property rather than a convenience. Variables not listed are unchanged.
    Bdd compose(const std::vector<std::pair<int, Bdd>> &subs) const;
    Bdd permute(const std::vector<std::pair<int, int>> &renaming) const;

    // --- structure ---------------------------------------------------------
    // top_var, low and high are meaningful only on a non-constant node and
    // throw otherwise. CUDD stores negation as a mark on the edge; low/high
    // resolve it, so a caller sees the plain Shannon decomposition.
    bool is_constant() const;
    bool is_true() const;
    bool is_false() const;
    int  top_var() const;
    Bdd  low() const;
    Bdd  high() const;

    // --- measurement -------------------------------------------------------
    // Satisfying assignments over a space of n_vars variables, as a double
    // because the honest answer routinely exceeds 64 bits.
    double      sat_count(int n_vars) const;
    std::size_t node_count() const;

private:
    Bdd(Manager *mgr, DdNode *node) noexcept : mgr_(mgr), node_(node) {}
    void release() noexcept;

    Manager *mgr_  = nullptr;
    DdNode  *node_ = nullptr;
};

} // namespace spbdd
