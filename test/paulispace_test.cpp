// ===========================================================================
//  PauliSpace builders, checked against the brute-force model in
//  reference.hpp.
//
//  The three builders that matter -- generated_by, weight_* and
//  commuting_with -- all reach their answer by a route the definition does not
//  suggest, so each is compared element by element against the definition on a
//  space small enough to enumerate.
// ===========================================================================

#include "check.hpp"
#include "reference.hpp"

#include <spbdd/spbdd.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace spbdd;

int main()
{
    SECTION("geometry");
    {
        PauliSpace sp(4);
        CHECK(sp.n_qubits() == 4 && sp.n_vars() == 8);
        CHECK(PauliSpace::xvar(3) == 6 && PauliSpace::zvar(3) == 7);
        CHECK_THROWS(PauliSpace(-1), std::invalid_argument);
    }

    SECTION("primitives");
    {
        PauliSpace sp(4);
        CHECK(sp.all().size() == 256.0);            // 4^4
        CHECK(sp.empty().is_empty());
        CHECK(sp.identity().size() == 1.0);
        CHECK(sp.identity().contains("IIII"));
        CHECK(sp.from("XYZI").size() == 1.0);
        CHECK(sp.from_list({"XXII", "IIZZ", "YYYY"}).size() == 3.0);
        CHECK(sp.parity({0, 2}, false).size() == 128.0);
        CHECK(sp.cube({{0, true}, {1, false}}).size() == 64.0);
        CHECK_THROWS(sp.from("XYZ"), std::invalid_argument);
        CHECK_THROWS(sp.from("XYZQ"), std::invalid_argument);
    }

    SECTION("support");
    {
        PauliSpace sp(4);
        CHECK(sp.supported_on({0, 1}).size() == 16.0);      // 4^2
        CHECK(sp.identity_on({0, 1}).size() == 16.0);
        CHECK(sp.pauli_at(2, Pauli::Y).size() == 64.0);     // 4^3
        CHECK(sp.pauli_at(2, 'Y') == sp.pauli_at(2, Pauli::Y));
        CHECK(sp.pauli_at(2, 'y') == sp.pauli_at(2, Pauli::Y));
        CHECK(sp.pauli_at(2, 'I') == sp.identity_on({2}));
        CHECK_THROWS(sp.pauli_at(2, 'Q'), std::invalid_argument);
        CHECK(sp.non_identity_at(2).size() == 192.0);       // 3 * 4^3
        CHECK(sp.matching("*X**").size() == 64.0);
        CHECK(sp.matching("IXYZ") == sp.from("IXYZ"));
        CHECK(sp.matching("****") == sp.all());
        CHECK_THROWS(sp.pauli_at(9, Pauli::X), std::out_of_range);
        CHECK_THROWS(sp.matching("**"), std::invalid_argument);
    }

    SECTION("weight vs the definition (n=4)");
    {
        PauliSpace sp(4);
        for (int w = 0; w <= 4; ++w)
            CHECK_AT(sp.weight_exactly(w) == spbdd_ref::as_set(sp, spbdd_ref::weight_between(4, w, w)),
                     ("weight_exactly(" + std::to_string(w) + ")").c_str());

        for (int w = 0; w <= 4; ++w)
            CHECK_AT(sp.weight_at_most(w) == spbdd_ref::as_set(sp, spbdd_ref::weight_between(4, 0, w)),
                     ("weight_at_most(" + std::to_string(w) + ")").c_str());

        CHECK(sp.weight_between(2, 3) == spbdd_ref::as_set(sp, spbdd_ref::weight_between(4, 2, 3)));

        // restricted to a subset of the qubits, the rest stay unconstrained
        const std::vector<int> qs = {0, 2};
        CHECK(sp.weight_exactly(1, qs) == spbdd_ref::as_set(sp, spbdd_ref::weight_between(4, 1, 1, qs)));
        CHECK(sp.weight_exactly(0, qs) == sp.identity_on(qs));
        CHECK(sp.weight_at_most(2, qs).is_universe());

        // out-of-range bounds are clamped rather than rejected
        CHECK(sp.weight_at_most(99) == sp.all());
        CHECK(sp.weight_between(5, 9).is_empty());
    }

    SECTION("generated_by vs the definition (n=3)");
    {
        PauliSpace sp(3);
        const std::vector<std::vector<std::string>> cases = {
            {},                              // no generators: the trivial group
            {"III"},                         // rank 0
            {"XXX", "XXX"},                  // dependent generators
            {"XXI", "IZZ"},
            {"XYZ", "ZZZ", "XXI"},
            {"XII", "IXI", "IIX", "ZII"},
        };
        for (const auto &gens : cases) {
            const std::string label = "generated_by(" + std::to_string(gens.size()) + " gens)";
            CHECK_AT(sp.generated_by(gens) == spbdd_ref::as_set(sp, spbdd_ref::span(gens, 3)),
                     label.c_str());
        }
        CHECK(sp.generated_by({}) == sp.identity());
    }

    SECTION("coset_of vs the definition (n=3)");
    {
        PauliSpace                     sp(3);
        const std::vector<std::string> gens = {"XXI", "IZZ"};
        for (const std::string base : {"III", "XYZ", "ZIX"})
            CHECK_AT(sp.coset_of(base, gens) == spbdd_ref::as_set(sp, spbdd_ref::coset(base, gens, 3)),
                     ("coset_of(" + base + ")").c_str());
        CHECK(sp.coset_of("III", gens) == sp.generated_by(gens));
    }

    SECTION("commutation vs the definition (n=3)");
    {
        PauliSpace sp(3);
        for (const std::string p : {"III", "XII", "XYZ", "ZZZ"}) {
            CHECK_AT(sp.commuting_with(p) == spbdd_ref::as_set(sp, spbdd_ref::commuting_with(p, 3)),
                     ("commuting_with(" + p + ")").c_str());
            CHECK_AT(sp.anticommuting_with(p) == ~sp.commuting_with(p),
                     ("anticommuting_with(" + p + ")").c_str());
        }
        CHECK(sp.commuting_with("III").is_universe());   // the identity commutes with everything

        const std::vector<std::string> gens = {"XXI", "IZZ"};
        CHECK(sp.commuting_with_all(gens) ==
              (sp.commuting_with("XXI") & sp.commuting_with("IZZ")));
    }

    SECTION("Steane code (n=7)");
    {
        PauliSpace                     sp(7);
        const std::vector<std::string> g = {"IIIXXXX", "IXXIIXX", "XIXIXIX",
                                            "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"};
        const PauliSet S = sp.generated_by(g);
        const PauliSet N = sp.commuting_with_all(g);

        CHECK(S.size() == 64.0);            // 2^rank, rank 6
        CHECK(N.size() == 256.0);           // 2^(14-6)
        CHECK(S.subset_of(N));              // the generators commute with each other
        CHECK((N - S).min_weight() == 3);   // the code distance
        CHECK((S - sp.identity()).min_weight() == 4);   // Steane's stabilisers are weight 4
    }

    SECTION("grow_to");
    {
        PauliSpace sp(2);
        PauliSpace copy = sp;               // handles share one universe
        const PauliSet f = sp.from("XZ");
        CHECK(f.size() == 1.0);

        sp.grow_to(4);
        CHECK(copy.n_qubits() == 4);        // growth is seen through every handle
        CHECK(f.size() == 16.0);            // now unconstrained on the new qubits
        CHECK((f & sp.identity_on({2, 3})).size() == 1.0);
        CHECK_THROWS(sp.grow_to(1), std::invalid_argument);
    }

    return REPORT("paulispace_test");
}
