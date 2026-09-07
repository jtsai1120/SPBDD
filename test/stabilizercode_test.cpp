// ===========================================================================
//  StabilizerCode: the symplectic basis and the inequivalent-pair decision.
//
//  The basis is checked against its defining relations, the syndrome and
//  signature maps against a second construction that never touches the
//  coordinate transform, and the decision procedure against brute force -- for
//  small codes, every pair of errors is actually formed and tested.
// ===========================================================================

#include "check.hpp"
#include "reference.hpp"

#include <spbdd/spbdd.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace spbdd;
using spbdd_ref::Strings;

namespace {

const std::vector<std::string> steane = {"IIIXXXX", "IXXIIXX", "XIXIXIX",
                                         "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"};
const std::vector<std::string> five_qubit = {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"};
const std::vector<std::string> bit_flip   = {"ZZI", "IZZ"};

bool commutes(const std::string &a, const std::string &b)
{
    return pauli_commute(parse_pauli_string(a), parse_pauli_string(b));
}

// The defining relations of Theorem "symplectic basis".
void check_basis(const StabilizerCode &code, const char *name)
{
    const auto &g = code.stabilizers();
    const auto &d = code.destabilizers();
    const auto &x = code.logical_x();
    const auto &z = code.logical_z();

    bool ok = true;
    for (std::size_t i = 0; i < g.size(); ++i) {
        for (std::size_t j = 0; j < g.size(); ++j) {
            if (!commutes(g[i], g[j])) ok = false;                       // <g,g> = 0
            if (!commutes(d[i], d[j])) ok = false;                       // <d,d> = 0
            if (commutes(g[i], d[j]) != (i != j)) ok = false;            // <g_i,d_j> = delta
        }
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (!commutes(g[i], x[j]) || !commutes(g[i], z[j])) ok = false;
            if (!commutes(d[i], x[j]) || !commutes(d[i], z[j])) ok = false;
        }
    }
    for (std::size_t i = 0; i < x.size(); ++i)
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (!commutes(x[i], x[j])) ok = false;
            if (!commutes(z[i], z[j])) ok = false;
            if (commutes(x[i], z[j]) != (i != j)) ok = false;             // <X,Z> = delta
        }

    CHECK_AT(ok, (std::string(name) + ": the symplectic basis relations hold").c_str());
}

// Brute force: form every pair and ask whether the product is in N(S) \ S.
bool brute_has_pair(const StabilizerCode &code, const Strings &errors)
{
    const PauliSet S = code.group();
    const PauliSet N = code.normalizer();
    for (const std::string &e1 : errors)
        for (const std::string &e2 : errors) {
            const std::string prod = pauli_string_to_text(
                pauli_mul(parse_pauli_string(e1), parse_pauli_string(e2)));
            if (N.contains(prod) && !S.contains(prod)) return true;
        }
    return false;
}

} // namespace

int main()
{
    SECTION("the symplectic basis");
    {
        PauliSpace s7(7), s5(5), s3(3);
        const StabilizerCode a(s7, steane);
        const StabilizerCode b(s5, five_qubit);
        const StabilizerCode c(s3, bit_flip);

        check_basis(a, "Steane [[7,1,3]]");
        check_basis(b, "five-qubit [[5,1,3]]");
        check_basis(c, "bit-flip [[3,1,1]]");

        CHECK(a.n_stabilizers() == 6 && a.n_logical() == 1);
        CHECK(b.n_stabilizers() == 4 && b.n_logical() == 1);
        CHECK(c.n_stabilizers() == 2 && c.n_logical() == 1);

        CHECK(a.group().size() == 64.0 && a.normalizer().size() == 256.0);
        CHECK(a.distance() == 3);
        CHECK(b.distance() == 3);
        CHECK(c.distance() == 1);   // the bit-flip code corrects no phase error

        // dependent generators are absorbed rather than rejected
        std::vector<std::string> repeated = steane;
        repeated.push_back(steane[0]);
        CHECK(StabilizerCode(s7, repeated).n_stabilizers() == 6);

        CHECK_THROWS(StabilizerCode(s3, {"XII", "ZII"}), std::invalid_argument);
        CHECK_THROWS(StabilizerCode(s3, {"XI"}), std::invalid_argument);
    }

    SECTION("syndrome and signature");
    {
        PauliSpace           sp(7);
        const StabilizerCode code(sp, steane);

        // The identity has the zero syndrome; the stabilizers do too, and they
        // carry no logical content either.
        CHECK(code.syndrome("IIIIIII") == std::vector<bool>(6, false));
        CHECK(code.logical_signature("IIIIIII") == std::vector<bool>(2, false));
        for (const std::string &g : code.stabilizers()) {
            CHECK_AT(code.syndrome(g) == std::vector<bool>(6, false), "a stabilizer has no syndrome");
            CHECK_AT(code.logical_signature(g) == std::vector<bool>(2, false),
                     "a stabilizer has no signature");
        }

        // Each destabilizer flips exactly its own syndrome bit; that is what
        // makes the pairing a symplectic basis.
        bool ok = true;
        for (int i = 0; i < 6; ++i) {
            const auto s = code.syndrome(code.destabilizers()[static_cast<std::size_t>(i)]);
            for (int j = 0; j < 6; ++j)
                if (s[static_cast<std::size_t>(j)] != (i == j)) ok = false;
        }
        CHECK(ok);

        // A logical operator is invisible to the syndrome but not to the
        // signature -- exactly the thing a decoder cannot see.
        const std::string lx = code.logical_x()[0];
        CHECK(code.syndrome(lx) == std::vector<bool>(6, false));
        CHECK(code.logical_signature(lx) != std::vector<bool>(2, false));
        CHECK(code.normalizer().contains(lx) && !code.group().contains(lx));

        // The class sets agree with the per-operator readings on every one of
        // the 4^7 operators, by construction of the sets rather than by
        // enumeration: each class is exactly the operators reading that value.
        std::vector<bool> s0(6, false), s1(6, false);
        s1[2] = true;
        CHECK(code.with_syndrome(s0).size() == 256.0);     // 4^7 / 2^6
        CHECK(code.with_syndrome(s1).size() == 256.0);
        CHECK(code.with_syndrome(s0).disjoint_from(code.with_syndrome(s1)));
        CHECK(code.with_syndrome(s0) == code.normalizer());   // zero syndrome is N(S)
        CHECK((code.with_syndrome(s0) & code.with_logical_signature({false, false})) ==
              code.group());                                  // ...and no logical part is S

        for (const std::string e : {"XIIIIII", "IZIIIII", "YYIIIII", "XXXXXXX"}) {
            CHECK_AT(code.with_syndrome(code.syndrome(e)).contains(e), "with_syndrome contains it");
            CHECK_AT(code.with_logical_signature(code.logical_signature(e)).contains(e),
                     "with_logical_signature contains it");
        }
    }

    SECTION("inequivalent pairs vs brute force (n=3)");
    {
        // Small enough to form every pair. The bit-flip code has distance 1, so
        // it is easy to find sets that fail.
        PauliSpace           sp(3);
        const StabilizerCode code(sp, bit_flip);

        const std::vector<Strings> cases = {
            {},
            {"III"},
            {"XII"},                        // one error is never a pair
            {"III", "ZZZ"},                 // product is a logical Z
            {"III", "ZII"},                 // product is a destabilizer: different syndrome
            {"XII", "IXI"},                 // same syndrome? different logical?
            {"III", "IIZ"},
            {"XII", "IXI", "IIX"},
            {"III", "XXX"},
        };

        for (const Strings &errors : cases) {
            const PauliSet    b     = spbdd_ref::as_set(sp, errors);
            const bool        got   = code.has_inequivalent_pair(b);
            const bool        want  = brute_has_pair(code, errors);
            const std::string label = "pair check on " + std::to_string(errors.size()) +
                                      " errors" + (want ? " (unsafe)" : " (safe)");
            CHECK_AT(got == want, label.c_str());
        }

        // Every subset of a five-element pool, checked exhaustively.
        const std::vector<std::string> pool = {"III", "XII", "IXI", "IIX", "ZZZ"};
        bool                           all_ok = true;
        for (int mask = 0; mask < 32; ++mask) {
            Strings errors;
            for (int i = 0; i < 5; ++i)
                if (mask & (1 << i)) errors.insert(pool[static_cast<std::size_t>(i)]);
            const PauliSet b = spbdd_ref::as_set(sp, errors);
            if (code.has_inequivalent_pair(b) != brute_has_pair(code, errors)) all_ok = false;
        }
        CHECK(all_ok);
    }

    SECTION("inequivalent pairs on the Steane code");
    {
        PauliSpace           sp(7);
        const StabilizerCode code(sp, steane);

        // Distance 3: any two errors of weight one are distinguishable, so a
        // set of them is safe.
        CHECK(!code.has_inequivalent_pair(sp.weight_at_most(1)));

        // Weight two reaches the distance, so it is not.
        CHECK(code.has_inequivalent_pair(sp.weight_at_most(2)));

        // The whole space certainly is, and a single error never is.
        CHECK(code.has_inequivalent_pair(sp.all()));
        CHECK(!code.has_inequivalent_pair(sp.from("XIIIIII")));
        CHECK(!code.has_inequivalent_pair(sp.empty()));

        // The stabilizer group itself is safe: every product stays inside S.
        CHECK(!code.has_inequivalent_pair(code.group()));

        // ...and N(S) is not, because a product can be a logical operator.
        CHECK(code.has_inequivalent_pair(code.normalizer()));
    }

    SECTION("witness extraction");
    {
        PauliSpace           sp(7);
        const StabilizerCode code(sp, steane);
        const PauliSet       errors = sp.weight_at_most(2);

        const auto found = code.find_inequivalent_pair(errors);
        CHECK(found.has_value());
        if (found) {
            const auto &p = *found;
            // Both witnesses are really in the set...
            CHECK(errors.contains(p.first));
            CHECK(errors.contains(p.second));
            // ...they really share a syndrome...
            CHECK(code.syndrome(p.first) == code.syndrome(p.second));
            CHECK(p.syndrome == code.syndrome(p.first));
            // ...their signatures really differ...
            CHECK(p.signature_first != p.signature_second);
            // ...and the product really is a logical operator.
            const std::string prod = pauli_string_to_text(
                pauli_mul(parse_pauli_string(p.first), parse_pauli_string(p.second)));
            CHECK(code.normalizer().contains(prod));
            CHECK(!code.group().contains(prod));

            std::printf("  ..    witnesses %s and %s, product %s\n",
                        p.first.c_str(), p.second.c_str(), prod.c_str());
        }

        CHECK(!code.find_inequivalent_pair(sp.weight_at_most(1)).has_value());
    }

    SECTION("the five-qubit code");
    {
        PauliSpace           sp(5);
        const StabilizerCode code(sp, five_qubit);

        CHECK(code.group().size() == 16.0);        // 2^4
        CHECK(code.normalizer().size() == 64.0);   // 2^(10-4)
        CHECK(code.distance() == 3);
        CHECK(!code.has_inequivalent_pair(sp.weight_at_most(1)));
        CHECK(code.has_inequivalent_pair(sp.weight_at_most(2)));

        const auto found = code.find_inequivalent_pair(sp.weight_at_most(2));
        CHECK(found.has_value());
        if (found)
            CHECK(code.syndrome(found->first) == code.syndrome(found->second) &&
                  found->signature_first != found->signature_second);
    }

    SECTION("argument checking and leaks");
    {
        PauliSpace           sp(7);
        const StabilizerCode code(sp, steane);
        PauliSpace           other(7);
        CHECK_THROWS(code.has_inequivalent_pair(other.all()), std::logic_error);
        CHECK_THROWS(code.with_syndrome({true}), std::invalid_argument);
        CHECK_THROWS(code.syndrome("XX"), std::invalid_argument);
    }

    {
        PauliSpace s5(5);
        {
            const StabilizerCode code(s5, five_qubit);
            (void)code.find_inequivalent_pair(s5.weight_at_most(2));
            (void)code.distance();
        }
        CHECK(s5.manager().check_zero_ref() == 0);
    }

    return REPORT("stabilizercode_test");
}
