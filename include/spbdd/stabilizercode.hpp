#pragma once

// ===========================================================================
//  spbdd::StabilizerCode -- the data a stabilizer code contributes to a query,
//  computed once and reused.
//
//  From a set of generators it builds a symplectic basis of GF(2)^(2n):
//
//      <g_i, g_j> = 0            <d_i, d_j> = 0            <g_i, d_j> = delta_ij
//      <g_i, Xbar_j> = <g_i, Zbar_j> = <d_i, Xbar_j> = <d_i, Zbar_j> = 0
//      <Xbar_i, Zbar_j> = delta_ij       <Xbar_i, Xbar_j> = <Zbar_i, Zbar_j> = 0
//
//  Every operator then expands as
//
//      e = sum a_i g_i + sum b_i d_i + sum (c_j Xbar_j + f_j Zbar_j)
//
//  and each coefficient is read by pairing with its *partner*, never with the
//  vector itself -- forced by the form being alternating:
//
//      a_i = <e, d_i>    b_i = <e, g_i>    c_j = <e, Zbar_j>   f_j = <e, Xbar_j>
//
//  b is the syndrome and (c, f) the logical signature. Two errors multiply
//  into a logical operator exactly when they share a syndrome and differ in
//  signature, which is what has_inequivalent_pair decides.
// ===========================================================================

#include "spbdd/pauliset.hpp"

#include <optional>
#include <string>
#include <vector>

namespace spbdd {

class StabilizerCode {
public:
    // The generators need not be independent, but they must commute with each
    // other; anything else is not a stabilizer group and throws.
    StabilizerCode(PauliSpace space, const std::vector<std::string> &generators);

    PauliSpace space() const { return sp_; }
    int        n_qubits() const { return sp_.n_qubits(); }
    int        n_stabilizers() const { return r_; }   // rank of the generator set
    int        n_logical() const { return k_; }       // n - r

    // --- the symplectic basis ----------------------------------------------
    const std::vector<std::string> &stabilizers() const { return g_; }
    const std::vector<std::string> &destabilizers() const { return d_; }
    const std::vector<std::string> &logical_x() const { return lx_; }
    const std::vector<std::string> &logical_z() const { return lz_; }

    // --- the sets ----------------------------------------------------------
    PauliSet group() const;         // S
    PauliSet normalizer() const;    // N(S)
    int      distance() const;      // the lightest element of N(S) \ S

    // --- reading one operator ----------------------------------------------
    // r bits, one per stabilizer: which of them the error anticommutes with.
    std::vector<bool> syndrome(const std::string &e) const;

    // 2k bits, (c_0..c_{k-1}, f_0..f_{k-1}).
    std::vector<bool> logical_signature(const std::string &e) const;

    // --- reading a whole class ---------------------------------------------
    // Built directly as parity checks in the original coordinates, so these do
    // not go through the transform below.
    PauliSet with_syndrome(const std::vector<bool> &syndrome) const;
    PauliSet with_logical_signature(const std::vector<bool> &signature) const;

    // --- the decision problem ----------------------------------------------
    // Does `errors` contain two operators whose product is a logical operator,
    // i.e. lies in N(S) \ S? Such a pair is indistinguishable to any decoder,
    // yet correcting one leaves a logical error in the other case.
    //
    // A constant number of diagram operations plus O(k), independent of how
    // many errors the set holds.
    bool has_inequivalent_pair(const PauliSet &errors) const;

    struct Pair {
        std::string       first, second;     // two members of `errors`
        std::vector<bool> syndrome;          // the one they share
        std::vector<bool> signature_first, signature_second;
    };

    // The same question, with a counterexample when the answer is yes.
    std::optional<Pair> find_inequivalent_pair(const PauliSet &errors) const;

private:
    using Row = std::vector<std::uint8_t>;   // a vector over GF(2), length 2n

    // Where each new coordinate lives among the 2n variables after the
    // transform. The blocks are contiguous: a, then b, then c, then f.
    int a_var(int i) const { return i; }
    int b_var(int i) const { return r_ + i; }
    int c_var(int j) const { return 2 * r_ + j; }
    int f_var(int j) const { return 2 * r_ + k_ + j; }

    // f . T^-1, i.e. the image of a set under the change of coordinates.
    // Points push forward but functions pull back, so it is the inverse that
    // gets substituted -- and T is not an involution, which is why it is
    // inverted explicitly at construction.
    Bdd to_code_coordinates(const Bdd &f) const;
    Row from_code_coordinates(const Row &y) const;   // e = T^-1 y

    PauliSpace sp_;
    int        r_ = 0;
    int        k_ = 0;

    std::vector<std::string> g_, d_, lx_, lz_;
    std::vector<Row>         grow_, drow_, lxrow_, lzrow_;   // the same, as vectors
    std::vector<Row>         tinv_;                          // T^-1, 2n x 2n
};

} // namespace spbdd
