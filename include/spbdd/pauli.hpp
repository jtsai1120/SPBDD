#pragma once

// ===========================================================================
//  Pauli operators in the symplectic (phase-free) representation.
//
//  A single-qubit Pauli is a pair of bits (x, z):
//
//      I = (0,0)   X = (1,0)   Z = (0,1)   Y = (1,1)
//
//  An n-qubit Pauli string is therefore a vector in GF(2)^(2n), and every
//  assignment of those 2n bits is a valid operator -- there is no unused
//  encoding space, which is why set complement needs no special handling
//  further up.
//
//  Phases are not tracked. That is a deliberate choice, not an omission:
//  without phases the group is abelian, multiplication is XOR, and the whole
//  of GF(2) linear algebra becomes available (it is what lets a subgroup be
//  described by parity checks instead of by listing its elements). The cost
//  is that P and -P are indistinguishable, which for questions like "does
//  this recovery fix that error" is the right identification anyway.
//
//  Nothing in this file knows what a BDD is.
// ===========================================================================

#include <cstdint>
#include <string>
#include <vector>

namespace spbdd {

// The enumerator's numeric value *is* the (x, z) bit pair, so no lookup table
// is ever needed between a Pauli and the two variables of its qubit.
enum class Pauli : std::uint8_t { I = 0, X = 1, Z = 2, Y = 3 };

constexpr bool has_x(Pauli p) { return (static_cast<std::uint8_t>(p) & 1u) != 0; }
constexpr bool has_z(Pauli p) { return (static_cast<std::uint8_t>(p) & 2u) != 0; }

constexpr Pauli make_pauli(bool x, bool z)
{
    return static_cast<Pauli>((x ? 1u : 0u) | (z ? 2u : 0u));
}

// 'I', 'X', 'Y', 'Z', either case. Anything else throws std::invalid_argument.
Pauli pauli_from_char(char c);
char  pauli_to_char(Pauli p);

// expected_n < 0 accepts any length; otherwise a mismatch throws
// std::invalid_argument, with `who` naming the caller in the message.
std::vector<Pauli> parse_pauli_string(const std::string &s, int expected_n = -1,
                                      const char *who = "parse_pauli_string");
std::string        pauli_string_to_text(const std::vector<Pauli> &p);

// Number of qubits on which the operator is not the identity.
int pauli_weight(const std::vector<Pauli> &p);

// True when the symplectic form of the two operators vanishes. Both must have
// the same length.
bool pauli_commute(const std::vector<Pauli> &a, const std::vector<Pauli> &b);

// Phase-free product, which over GF(2) is exactly XOR.
Pauli              pauli_mul(Pauli a, Pauli b);
std::vector<Pauli> pauli_mul(const std::vector<Pauli> &a, const std::vector<Pauli> &b);

} // namespace spbdd
