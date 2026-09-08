#pragma once

// ===========================================================================
//  spbdd::Bdd -- an owning handle on one BDD node.
//
//  BuDDy's own header is deliberately absent here. A BuDDy node is an index
//  into a global table rather than a pointer, so the handle below is an int
//  and this header needs nothing from <bdd.h>: a program that uses spbdd needs
//  -Iinclude and nothing else on its include path. The price is that every
//  function which has to call into BuDDy is defined in src/bdd.cpp instead of
//  inline; next to the cost of a BDD operation that is not measurable.
//
//  Two properties everything above this layer relies on:
//
//    1. Reference counting. A Bdd owns exactly one reference to its node and
//       hands it back when it dies. No code above this file writes a
//       bdd_addref by hand.
//
//    2. Canonicity. Two Bdds are equal iff they denote the same boolean
//       function, so equality is an integer comparison, never a traversal.
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

namespace spbdd {

class Manager;

// BuDDy's node handle. Plain `int`, so no forward declaration is possible or
// needed; the zero value is a legitimate node (the false terminal), which is
// why validity is carried by the Manager pointer rather than by this.
using NodeId = int;

class Bdd {
public:
    Bdd() noexcept = default;

    // Takes over a reference the caller already owns -- that is, exactly what
    // a BuDDy call returns once it has been bdd_addref'd.
    static Bdd adopt(Manager *mgr, NodeId node) noexcept;

    // Adds a reference of its own.
    static Bdd borrow(Manager *mgr, NodeId node);

    Bdd(const Bdd &);
    Bdd(Bdd &&) noexcept;
    Bdd &operator=(const Bdd &);
    Bdd &operator=(Bdd &&) noexcept;
    ~Bdd();

    bool     valid() const noexcept { return mgr_ != nullptr; }
    Manager *manager() const noexcept { return mgr_; }
    NodeId   node() const noexcept { return node_; }

    // Canonicity (see the header comment) makes this an integer comparison.
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
    // `vars` is a plain list of variable numbers; the varset BuDDy wants is
    // built and released internally.
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
    // throw otherwise. BuDDy has no complement edges, so these are already the
    // plain Shannon decomposition with nothing to resolve.
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

    // Nodes in the diagram, counting the terminal -- so a constant is 1 and
    // "x0 and x1" is 3. BuDDy counts internal nodes only, so this adds one.
    std::size_t node_count() const;

private:
    Bdd(Manager *mgr, NodeId node) noexcept : mgr_(mgr), node_(node) {}
    void release() noexcept;

    Manager *mgr_  = nullptr;
    NodeId   node_ = 0;
};

} // namespace spbdd
