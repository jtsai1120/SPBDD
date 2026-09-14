#pragma once

// ===========================================================================
//  spbdd::PauliSet -- one set of Pauli operators inside a PauliSpace.
//
//  Two members: the space it lives in, and the characteristic function of the
//  set as a BDD. The space is held by value because PauliSpace *is* a handle
//  (one shared_ptr), so holding it costs a pointer and keeps the manager alive
//  for exactly as long as some set still refers to it.
//
//  Operands of a binary operation must come from the same space.
// ===========================================================================

#include "spbdd/paulispace.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace spbdd {

// The two halves a measurement splits a set into; defined below PauliSet,
// which it holds by value.
struct MeasurementSplit;

class PauliSet {
public:
    PauliSpace space() const { return sp_; }
    int        n_qubits() const { return sp_.n_qubits(); }

    const Bdd &bdd() const { return f_; }        // escape hatch

    // --- set algebra -------------------------------------------------------
    // Every one of the 4^n assignments is a valid operator, so the complement
    // is the plain boolean negation with no universe to intersect against.
    PauliSet operator|(const PauliSet &o) const;   // union
    PauliSet operator&(const PauliSet &o) const;   // intersection
    PauliSet operator-(const PauliSet &o) const;   // difference
    PauliSet operator^(const PauliSet &o) const;   // symmetric difference
    PauliSet operator~() const;                    // complement

    // --- product -----------------------------------------------------------
    // {a b : a in this, b in o}: every element of one set multiplied by every
    // element of the other. Without phases the product is the XOR of the two
    // symplectic vectors, so this is commutative, and it is what composing
    // two sources of error means: fault_inject(Q) is the product with
    // supported_on(Q), a coset is the product with a group, and two elements
    // of E multiply into a logical operator iff E * E meets N(S) \ S. The cost
    // is bounded by the product of the two diagram sizes, never by the number
    // of elements.
    PauliSet  operator*(const PauliSet &o) const;
    PauliSet &operator*=(const PauliSet &o);

    PauliSet &operator|=(const PauliSet &o);
    PauliSet &operator&=(const PauliSet &o);
    PauliSet &operator-=(const PauliSet &o);
    PauliSet &operator^=(const PauliSet &o);

    // Canonicity makes this a node comparison rather than a traversal.
    bool operator==(const PauliSet &o) const;
    bool operator!=(const PauliSet &o) const { return !(*this == o); }

    // --- predicates --------------------------------------------------------
    bool is_empty() const;
    bool is_universe() const;
    bool contains(const std::string &p) const;
    bool contains(const std::vector<Pauli> &p) const;
    bool subset_of(const PauliSet &o) const;
    bool superset_of(const PauliSet &o) const;
    bool disjoint_from(const PauliSet &o) const;

    // --- metrics -----------------------------------------------------------
    // The cardinality, which reaches 4^n and so is a double. Unrelated to
    // node_count(), which is what the set actually costs.
    double      size() const;
    std::size_t node_count() const;

    // The smallest weight present, or -1 for the empty set. Binary search over
    // the weight bound, so log(n) diagram constructions and no enumeration.
    int min_weight() const;

    // --- derived sets ------------------------------------------------------
    // Multiply by every Pauli supported on these qubits. That is one
    // existential quantification, not 4^k multiplications: the operators being
    // applied are exactly the coordinate subspace spanned by those axes.
    PauliSet forget(const std::vector<int> &qubits) const;

    // The same operation under the name the circuit level uses: a fault at
    // these locations can leave behind any Pauli supported on them.
    PauliSet fault_inject(const std::vector<int> &qubits) const;

    // Rewrite every element to be the identity on these qubits. Quantify
    // first, then pin -- the other order deletes information that the
    // quantification would have recovered.
    PauliSet reset(const std::vector<int> &qubits) const;

    PauliSet with_weight_at_most(int w) const;
    PauliSet with_weight_exactly(int w) const;

    // --- commutation -------------------------------------------------------
    // The elements that commute, or anticommute, with a Pauli operator. The
    // two are complementary within the set, and each is one intersection with
    // a half-space. Measuring an observable `p` on a set of errors is exactly
    // this split: the commuting half sees the outcome the ideal circuit would
    // give, the anticommuting half sees it flipped, so measuring a stabilizer
    // generator is one bit of syndrome extraction.
    PauliSet with_commuting(const std::string &p) const;
    PauliSet with_anticommuting(const std::string &p) const;

    // --- measurement -------------------------------------------------------
    // The one-qubit cases, returning both halves at once. An empty side is an
    // outcome that cannot occur under these errors. The measurement does not
    // by itself do anything to the qubit; a destructive measurement that reuses
    // the qubit is this split followed by reset().
    MeasurementSplit measure_z(int qubit) const;   // reads the x coordinate
    MeasurementSplit measure_x(int qubit) const;   // reads the z coordinate
    MeasurementSplit measure_y(int qubit) const;   // reads their parity

    // --- Clifford gates ----------------------------------------------------
    // Every gate acts by conjugation, E -> U E U*. The elements of a set are
    // deviations from an ideal state, and applying U moves the ideal state
    // too, so relative to the new ideal state the residual error is U E U*.
    // There is no second semantics.
    //
    // Each of these is a linear map on GF(2)^(2n), applied as one simultaneous
    // substitution -- a map like CX's {x_t := x_t ^ x_c, z_c := z_c ^ z_t}
    // read one assignment at a time would give the wrong answer.
    //
    // All of them are involutions here (their squares are Pauli operators,
    // which act trivially without phases), so substituting the map and
    // substituting its inverse coincide. A gate whose square is not a Pauli
    // would have to substitute the inverse.

    // Conjugation by a Pauli only ever changes a sign, so these three are the
    // identity map. They exist so a circuit can be replayed without having to
    // special-case them.
    PauliSet x(int qubit) const;
    PauliSet y(int qubit) const;
    PauliSet z(int qubit) const;

    PauliSet h(int qubit) const;    // X <-> Z
    PauliSet s(int qubit) const;    // X -> Y, Z -> Z.  S and S* agree here
    PauliSet sdg(int qubit) const;  // S*: same map as s() in this representation
    PauliSet sx(int qubit) const;   // sqrt(X):  X -> X, Z -> Y
    PauliSet sxdg(int qubit) const; // sqrt(X)*: same map as sx() in this representation

    PauliSet cx(int control, int target) const;
    PauliSet cy(int control, int target) const;
    PauliSet cz(int a, int b) const;
    PauliSet swap(int a, int b) const;

    // --- enumeration -------------------------------------------------------
    // These are the operations whose cost really is the number of elements.
    //
    // for_each calls fn once per element, and stops early if fn returns false.
    // The order depends on the current variable order and is not specified.
    void for_each(const std::function<bool(const std::string &)> &fn) const;
    void for_each(const std::function<bool(const std::vector<Pauli> &)> &fn) const;

    // Sorted, so the result does not depend on the variable order. Throws
    // std::length_error rather than quietly producing a huge vector.
    std::vector<std::string> to_strings(std::size_t max_elements = 4096) const;

    // Any one element, found by descending the diagram once: O(n), whatever
    // the cardinality.
    std::optional<std::string> any_element() const;

private:
    friend class PauliSpace;
    PauliSet(PauliSpace sp, Bdd f) : sp_(std::move(sp)), f_(std::move(f)) {}

    PauliSpace sp_;
    Bdd        f_;
};

struct MeasurementSplit {
    PauliSet unflipped;   // the outcome the ideal circuit would give
    PauliSet flipped;     // the other one
};

} // namespace spbdd
