#pragma once

// ===========================================================================
//  spbdd::PauliSpace -- the universe of n-qubit Pauli operators, and the
//  factory for every set built inside it.
//
//  A set of Pauli operators is stored as the characteristic function of the
//  corresponding subset of GF(2)^(2n), as a BDD over 2n boolean variables.
//  The consequence that makes the whole library worthwhile is that the size of
//  the diagram and the cardinality of the set are unrelated: a set of 4^n
//  elements can have a handful of nodes. So no builder here is allowed to work
//  by enumerating elements, and none does.
//
//  The variable convention lives in this file and nowhere else:
//
//      xvar(q) = 2q        zvar(q) = 2q + 1
//
//  It is only the *declaration* order. CUDD may reorder levels freely; all
//  code here is written against variable numbers, which never change.
//
//  Nothing in this file knows that CUDD exists.
//
//  Note that every builder returns a PauliSet, which is declared in
//  pauliset.hpp -- so this header compiles on its own but cannot be used on
//  its own. Include spbdd.hpp.
// ===========================================================================

#include "spbdd/manager.hpp"
#include "spbdd/pauli.hpp"

#include <initializer_list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace spbdd {

class PauliSet;

// A copyable handle on one universe. Copies share the same manager, which
// stays alive until the last PauliSpace and PauliSet referring to it is gone.
class PauliSpace {
public:
    explicit PauliSpace(int n_qubits, const ManagerConfig &cfg = {});

    int n_qubits() const;
    int n_vars() const;                          // 2 * n_qubits

    static int xvar(int q) { return 2 * q; }
    static int zvar(int q) { return 2 * q + 1; }

    Manager &manager() const;                    // escape hatch

    // Only grows. Sets built before the growth are unchanged as functions,
    // which means they now place no constraint on the new qubits -- so an
    // element of such a set becomes 4^k elements, one per Pauli on the k new
    // qubits. Intersect with identity_on(new qubits) if that is not wanted.
    void grow_to(int n_qubits);

    // --- primitives --------------------------------------------------------
    PauliSet empty() const;
    PauliSet all() const;
    PauliSet identity() const;
    PauliSet literal(int var, bool positive) const;
    PauliSet cube(const std::vector<std::pair<int, bool>> &assignment) const;

    // Everything satisfying "the XOR of these variables equals target". This
    // is the workhorse: subgroups and commutation conditions are both just
    // conjunctions of these, and each one costs about 2n nodes.
    PauliSet parity(const std::vector<int> &vars, bool target) const;

    // Wrap a Bdd built by hand through manager().
    PauliSet wrap(Bdd f) const;

    // --- explicitly listed elements ----------------------------------------
    PauliSet from(const std::string &s) const;
    PauliSet from(const std::vector<Pauli> &p) const;
    // The union of the listed operators.
    PauliSet from(const std::vector<std::string> &list) const;
    // A braced list of string literals also matches the iterator-pair
    // constructors of std::string and std::vector<Pauli>, which would make
    // from({XXII, IIZZ}) ambiguous. This exact match settles it.
    PauliSet from(std::initializer_list<std::string> list) const;

    // --- subgroups ---------------------------------------------------------
    // The GF(2) span of the generators: 2^rank elements, built without ever
    // forming a product. See the note in paulispace.cpp.
    PauliSet generated_by(const std::vector<std::string> &generators) const;
    PauliSet coset_of(const std::string &base, const std::vector<std::string> &generators) const;

    // --- support -----------------------------------------------------------
    PauliSet supported_on(const std::vector<int> &qubits) const;   // identity elsewhere
    PauliSet identity_on(const std::vector<int> &qubits) const;    // free elsewhere
    PauliSet pauli_at(int qubit, Pauli p) const;
    PauliSet pauli_at(int qubit, char p) const;   // 'I' 'X' 'Y' 'Z', either case
    PauliSet non_identity_at(int qubit) const;

    // Glob over the n qubits: I X Y Z pin a qubit, '*' '?' '.' leave it free.
    PauliSet matching(const std::string &pattern) const;

    // --- weight ------------------------------------------------------------
    // Built with a running counter, so the cost is O(n*w) nodes rather than
    // the C(n,w) cubes the definition suggests.
    PauliSet weight_exactly(int w) const;
    PauliSet weight_at_most(int w) const;
    PauliSet weight_between(int lo, int hi) const;

    // The same, counting only the listed qubits and leaving the rest free.
    PauliSet weight_exactly(int w, const std::vector<int> &qubits) const;
    PauliSet weight_at_most(int w, const std::vector<int> &qubits) const;
    PauliSet weight_between(int lo, int hi, const std::vector<int> &qubits) const;

    // --- commutation -------------------------------------------------------
    PauliSet commuting_with(const std::string &p) const;
    PauliSet anticommuting_with(const std::string &p) const;

    // Everything commuting with all of them at once. For the generators of a
    // stabilizer group S this is the normaliser N(S).
    PauliSet commuting_with_all(const std::vector<std::string> &ps) const;

private:
    // n_qubits lives beside the manager rather than in the handle, so that
    // grow_to() is seen by every copy.
    struct State {
        Manager mgr;
        int     n;
        State(int n_qubits, const ManagerConfig &cfg) : mgr(2 * n_qubits, cfg), n(n_qubits) {}
    };

    std::shared_ptr<State> st_;
};

} // namespace spbdd
