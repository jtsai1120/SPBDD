// ===========================================================================
//  PauliSet: set algebra, predicates, metrics, derived sets and enumeration,
//  checked against the brute-force model in reference.hpp.
// ===========================================================================

#include "check.hpp"
#include "reference.hpp"

#include <spbdd/spbdd.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

using namespace spbdd;
using spbdd_ref::Strings;

static void checks_on_three_qubits()
{
    // Two arbitrary but overlapping subsets of the 64 operators on 3 qubits.
    PauliSpace sp(3);
    Strings    A, B;
    {
        const auto all = spbdd_ref::universe(3);
        int        i   = 0;
        for (const std::string &s : all) {
            if (i % 3 == 0) A.insert(s);
            if (i % 5 == 0 || i % 3 == 0) B.insert(s);
            ++i;
        }
    }
    const PauliSet a = spbdd_ref::as_set(sp, A);
    const PauliSet b = spbdd_ref::as_set(sp, B);

    SECTION("set algebra vs the definition (n=3)");
    {
        CHECK((a | b) == spbdd_ref::as_set(sp, spbdd_ref::unite(A, B)));
        CHECK((a & b) == spbdd_ref::as_set(sp, spbdd_ref::intersect(A, B)));
        CHECK((a - b) == spbdd_ref::as_set(sp, spbdd_ref::subtract(A, B)));
        CHECK((a ^ b) == spbdd_ref::as_set(sp, spbdd_ref::sym_diff(A, B)));
        CHECK((~a) == spbdd_ref::as_set(sp, spbdd_ref::complement(A, 3)));
        CHECK((a | ~a).is_universe());
        CHECK((a & ~a).is_empty());

        // the compound forms agree with the plain ones
        PauliSet acc = a;
        acc |= b;
        acc &= b;
        acc -= a;
        acc ^= b;
        CHECK(acc == ((((a | b) & b) - a) ^ b));

        // BuDDy keeps its tables in process-wide globals, so a second space
        // cannot be built at all. That takes the place of the "operands from
        // different spaces" check, which has nothing to run against here.
        CHECK_THROWS(PauliSpace(3), std::runtime_error);
    }

    SECTION("predicates");
    {
        CHECK(sp.empty().is_empty() && !sp.empty().is_universe());
        CHECK(sp.all().is_universe() && !sp.all().is_empty());
        CHECK(a.subset_of(b) && b.superset_of(a) && !b.subset_of(a));
        CHECK(a.subset_of(a));
        CHECK(a.disjoint_from(~a) && !a.disjoint_from(a));
        CHECK(sp.empty().disjoint_from(a));
        CHECK(a.contains(*A.begin()));
        CHECK(!a.contains(*spbdd_ref::complement(A, 3).begin()));
        CHECK_THROWS(a.contains("XX"), std::invalid_argument);
    }

    SECTION("metrics");
    {
        CHECK(a.size() == static_cast<double>(A.size()));
        CHECK(sp.all().size() == 64.0);
        CHECK(sp.empty().size() == 0.0);
        CHECK(sp.empty().node_count() == 1);      // the constant node
        CHECK(a.node_count() > 1);
    }

    SECTION("min_weight vs enumeration");
    {
        // against the definition: the lightest element found by enumeration
        for (const std::string p : {"IXI", "ZZZ", "XYZ"}) {
            const PauliSet s    = sp.commuting_with(p) - sp.identity();
            int            want = 4;
            s.for_each([&](const std::vector<Pauli> &e) {
                want = std::min(want, pauli_weight(e));
                return true;
            });
            CHECK_AT(s.min_weight() == want, ("min_weight vs enumeration, " + p).c_str());
        }
    }

    SECTION("forget / reset vs the definition (n=3)");
    {
        const Strings          elems = {"XIZ", "IYI", "ZZZ"};
        const PauliSet         set   = spbdd_ref::as_set(sp, elems);
        const std::vector<int> Q     = {1};

        CHECK(set.forget(Q) == spbdd_ref::as_set(sp, spbdd_ref::forget(elems, Q, 3)));
        CHECK(set.reset(Q) == spbdd_ref::as_set(sp, spbdd_ref::reset(elems, Q)));
        CHECK(set.forget(Q).size() == 12.0);      // 3 elements * 4 Paulis, all distinct
        CHECK(set.reset(Q).size() == 3.0);
        CHECK(set.subset_of(set.forget(Q)));      // quantification only ever adds
        CHECK(set.forget({}) == set);
        CHECK(set.fault_inject(Q) == set.forget(Q));   // the circuit-level name
        CHECK(set.fault_inject({0, 2}) == set.forget({0, 2}));
        CHECK(set.reset(Q).subset_of(sp.identity_on(Q)));
        CHECK(sp.all().forget({0, 1, 2}) == sp.all());
        CHECK(sp.all().reset({0, 1, 2}) == sp.identity());
        CHECK_THROWS(set.forget({9}), std::out_of_range);

        // derived weight sets are just intersections
        CHECK(a.with_weight_at_most(1) == (a & sp.weight_at_most(1)));
        CHECK(a.with_weight_exactly(2) == (a & sp.weight_exactly(2)));
    }

    SECTION("measurement");
    {
        // The two sides partition the set: an element either commutes with the
        // observable or it does not.
        for (const std::string obs : {"ZII", "IXI", "YYY", "XZI", "III"}) {
            const MeasurementSplit m = a.measure(obs);
            const bool             ok =
                (m.unflipped | m.flipped) == a && m.unflipped.disjoint_from(m.flipped);
            CHECK_AT(ok, ("measure(" + obs + ") partitions the set").c_str());

            // and each side is exactly what the definition says it is
            Strings want_unflipped, want_flipped;
            for (const std::string &e : A) {
                if (pauli_commute(parse_pauli_string(e), parse_pauli_string(obs)))
                    want_unflipped.insert(e);
                else
                    want_flipped.insert(e);
            }
            const bool matches = m.unflipped == spbdd_ref::as_set(sp, want_unflipped) &&
                                 m.flipped == spbdd_ref::as_set(sp, want_flipped);
            CHECK_AT(matches, ("measure(" + obs + ") matches the definition").c_str());
        }

        // Measuring the identity can never flip anything.
        CHECK(a.measure("III").flipped.is_empty());
        CHECK(a.measure("III").unflipped == a);

        // A Z-basis measurement reads the x coordinate, an X-basis measurement
        // the z coordinate, and a Y-basis one their parity.
        CHECK(sp.all().measure_z(0).flipped == sp.all().measure("ZII").flipped);
        CHECK(sp.all().measure_x(1).flipped == sp.all().measure("IXI").flipped);
        CHECK(sp.all().measure_y(2).flipped == sp.all().measure("IIY").flipped);
        CHECK(sp.from("XII").measure_z(0).flipped.size() == 1.0);
        CHECK(sp.from("ZII").measure_z(0).flipped.is_empty());
        CHECK(sp.from("YII").measure_z(0).flipped.size() == 1.0);
        CHECK(sp.from("XII").measure_x(0).flipped.is_empty());
        CHECK(sp.from("YII").measure_y(0).flipped.is_empty());

        // Half the universe on each side of a non-trivial observable.
        CHECK(sp.all().measure_z(0).unflipped.size() == 32.0);
        CHECK(sp.all().measure_z(0).flipped.size() == 32.0);

        CHECK_THROWS(sp.all().measure_z(9), std::out_of_range);
        CHECK_THROWS(sp.all().measure("XX"), std::invalid_argument);
    }
}

// Six qubits, where min_weight has room to bisect. A separate function because
// only one PauliSpace may be alive at a time on this backend.
static void checks_on_six_qubits()
{
    SECTION("min_weight");
    {
        PauliSpace s6(6);
        CHECK(s6.empty().min_weight() == -1);
        CHECK(s6.identity().min_weight() == 0);
        CHECK(s6.all().min_weight() == 0);
        CHECK(s6.from("XXIIII").min_weight() == 2);
        CHECK(s6.from_list({"XXXXXI", "IIZZII"}).min_weight() == 2);
        CHECK(s6.weight_exactly(4).min_weight() == 4);
        CHECK(s6.weight_exactly(6).min_weight() == 6);
    }
}

int main()
{
    checks_on_three_qubits();
    checks_on_six_qubits();

    SECTION("enumeration");
    {
        PauliSpace     s4(4);
        const PauliSet w = s4.weight_exactly(2);

        // the walk and Cudd_CountMinterm are independent routes to the count
        std::size_t seen = 0;
        w.for_each([&](const std::vector<Pauli> &) { ++seen; return true; });
        CHECK(static_cast<double>(seen) == w.size());

        // the string form walks the same elements and stops the same way
        std::vector<std::string> via_string;
        w.for_each([&](const std::string &e) { via_string.push_back(e); return true; });
        std::sort(via_string.begin(), via_string.end());
        CHECK(via_string == w.to_strings());

        std::size_t string_stopped = 0;
        w.for_each([&](const std::string &) { return ++string_stopped < 3; });
        CHECK(string_stopped == 3);

        std::size_t stopped = 0;
        w.for_each([&](const std::vector<Pauli> &) { return ++stopped < 5; });
        CHECK(stopped == 5);

        const auto got = w.to_strings();
        CHECK(got == spbdd_ref::as_vector(spbdd_ref::weight_between(4, 2, 2)));
        CHECK(std::is_sorted(got.begin(), got.end()));

        // sorting is what makes this independent of the variable order
        s4.manager().reorder_now();
        CHECK(w.to_strings() == got);

        CHECK_THROWS(w.to_strings(3), std::length_error);
        CHECK(s4.empty().to_strings().empty());
        CHECK(s4.identity().to_strings() == std::vector<std::string>{"IIII"});
        CHECK(s4.all().to_strings(1000).size() == 256u);

        CHECK(!s4.empty().any_element().has_value());
        for (const PauliSet &s : {w, s4.all(), s4.identity(), s4.from("XYZI")}) {
            const auto e = s.any_element();
            CHECK_AT(e.has_value() && s.contains(*e), "any_element is a member");
        }
    }

    SECTION("reference counting");
    {
        PauliSpace s5(5);
        {
            const PauliSet x = s5.weight_at_most(2);
            const PauliSet y = s5.commuting_with("XZXZX");
            const PauliSet z = (x & y) | (~x ^ y);
            (void)z.to_strings(100000);
            (void)z.min_weight();
            (void)z.forget({0, 3});
        }
        CHECK(s5.manager().check_zero_ref() == 0);
    }

    return REPORT("pauliset_test");
}
