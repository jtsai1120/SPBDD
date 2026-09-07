#include "spbdd/pauli.hpp"

#include <stdexcept>
#include <string>

namespace spbdd {

Pauli pauli_from_char(char c)
{
    switch (c) {
    case 'I': case 'i': return Pauli::I;
    case 'X': case 'x': return Pauli::X;
    case 'Y': case 'y': return Pauli::Y;
    case 'Z': case 'z': return Pauli::Z;
    default:
        throw std::invalid_argument(std::string("not a Pauli letter: '") + c + "'");
    }
}

char pauli_to_char(Pauli p)
{
    switch (p) {
    case Pauli::I: return 'I';
    case Pauli::X: return 'X';
    case Pauli::Z: return 'Z';
    case Pauli::Y: return 'Y';
    }
    return '?';   // unreachable while Pauli holds one of its four enumerators
}

std::vector<Pauli> parse_pauli_string(const std::string &s, int expected_n, const char *who)
{
    if (expected_n >= 0 && static_cast<int>(s.size()) != expected_n)
        throw std::invalid_argument(std::string(who) + ": expected " +
                                    std::to_string(expected_n) + " qubits, got \"" + s +
                                    "\" of length " + std::to_string(s.size()));

    std::vector<Pauli> out;
    out.reserve(s.size());
    for (char c : s) out.push_back(pauli_from_char(c));
    return out;
}

std::string pauli_string_to_text(const std::vector<Pauli> &p)
{
    std::string s;
    s.reserve(p.size());
    for (Pauli q : p) s.push_back(pauli_to_char(q));
    return s;
}

int pauli_weight(const std::vector<Pauli> &p)
{
    int w = 0;
    for (Pauli q : p)
        if (q != Pauli::I) ++w;
    return w;
}

bool pauli_commute(const std::vector<Pauli> &a, const std::vector<Pauli> &b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("pauli_commute: operands have different lengths");

    // The symplectic form: sum over qubits of x(a)z(b) + z(a)x(b), mod 2.
    unsigned parity = 0;
    for (std::size_t q = 0; q < a.size(); ++q)
        parity ^= static_cast<unsigned>(has_x(a[q]) && has_z(b[q])) ^
                  static_cast<unsigned>(has_z(a[q]) && has_x(b[q]));
    return parity == 0;
}

Pauli pauli_mul(Pauli a, Pauli b)
{
    return static_cast<Pauli>(static_cast<std::uint8_t>(a) ^ static_cast<std::uint8_t>(b));
}

std::vector<Pauli> pauli_mul(const std::vector<Pauli> &a, const std::vector<Pauli> &b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("pauli_mul: operands have different lengths");

    std::vector<Pauli> out;
    out.reserve(a.size());
    for (std::size_t q = 0; q < a.size(); ++q) out.push_back(pauli_mul(a[q], b[q]));
    return out;
}

} // namespace spbdd
