// ===========================================================================
//  Clifford gates.
//
//  Each gate is checked three ways:
//
//    * against the generator table, applied element by element by the
//      reference model -- the physics, on a space small enough to enumerate;
//    * as a map: it must be a bijection, and (here) an involution;
//    * as a symplectic map: it must preserve commutation, which is what makes
//      the image of a stabilizer group a stabilizer group.
// ===========================================================================

#include "check.hpp"
#include "reference.hpp"

#include <spbdd/spbdd.hpp>

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace spbdd;
using spbdd_ref::Strings;

namespace {

// The generator table for every gate, in the order X, Z of each named qubit.
struct GateSpec {
    const char                     *name;
    std::vector<int>                qubits;
    std::vector<std::string>        images;
    std::function<PauliSet(const PauliSet &)> apply;
};

} // namespace

int main()
{
    // Only one PauliSpace may be alive at a time on this backend, so each
    // section below builds its own rather than sharing one across main.
    const int n = 3;

    // Qubit 2 is deliberately left out of every gate, so the tests also see
    // that untouched qubits are untouched.
    const std::vector<GateSpec> gates = {
        {"x(0)",     {0},    {"X", "Z"},                    [](const PauliSet &s) { return s.x(0); }},
        {"y(0)",     {0},    {"X", "Z"},                    [](const PauliSet &s) { return s.y(0); }},
        {"z(0)",     {0},    {"X", "Z"},                    [](const PauliSet &s) { return s.z(0); }},
        {"h(0)",     {0},    {"Z", "X"},                    [](const PauliSet &s) { return s.h(0); }},
        {"h(1)",     {1},    {"Z", "X"},                    [](const PauliSet &s) { return s.h(1); }},
        {"s(0)",     {0},    {"Y", "Z"},                    [](const PauliSet &s) { return s.s(0); }},
        {"sx(0)",    {0},    {"X", "Y"},                    [](const PauliSet &s) { return s.sx(0); }},
        {"cx(0,1)",  {0, 1}, {"XX", "ZI", "IX", "ZZ"},      [](const PauliSet &s) { return s.cx(0, 1); }},
        {"cx(1,0)",  {1, 0}, {"XX", "ZI", "IX", "ZZ"},      [](const PauliSet &s) { return s.cx(1, 0); }},
        {"cy(0,1)",  {0, 1}, {"XY", "ZI", "ZX", "ZZ"},      [](const PauliSet &s) { return s.cy(0, 1); }},
        {"cz(0,1)",  {0, 1}, {"XZ", "ZI", "ZX", "IZ"},      [](const PauliSet &s) { return s.cz(0, 1); }},
        {"swap(0,1)",{0, 1}, {"IX", "IZ", "XI", "ZI"},      [](const PauliSet &s) { return s.swap(0, 1); }},
    };

    SECTION("generator table, element by element (n=3)");
    {
        PauliSpace sp(n);
        // The whole universe at once: 64 operators, every one of them checked.
        const Strings  all     = spbdd_ref::universe(n);
        const PauliSet all_set = sp.all();

        for (const GateSpec &g : gates) {
            const Strings want = spbdd_ref::conjugate(all, g.qubits, g.images, n);
            CHECK_AT(g.apply(all_set).is_universe(), (std::string(g.name) + " maps the universe onto itself").c_str());
            CHECK_AT(want.size() == all.size(), (std::string(g.name) + " is injective").c_str());

            // and on a handful of individual operators
            bool ok = true;
            for (const std::string &e : all) {
                const Strings one  = {e};
                const Strings img  = spbdd_ref::conjugate(one, g.qubits, g.images, n);
                if (g.apply(sp.from(e)) != spbdd_ref::as_set(sp, img)) { ok = false; break; }
            }
            CHECK_AT(ok, (std::string(g.name) + " matches the generator table on all 64").c_str());
        }
    }

    SECTION("named images");
    {
        PauliSpace sp(n);
        // The table above spelled out; these read as the physics.
        CHECK(sp.from("XII").h(0) == sp.from("ZII"));
        CHECK(sp.from("ZII").h(0) == sp.from("XII"));
        CHECK(sp.from("YII").h(0) == sp.from("YII"));
        CHECK(sp.from("XII").s(0) == sp.from("YII"));
        CHECK(sp.from("ZII").s(0) == sp.from("ZII"));
        CHECK(sp.from("ZII").sx(0) == sp.from("YII"));
        CHECK(sp.from("XII").cx(0, 1) == sp.from("XXI"));
        CHECK(sp.from("IXI").cx(0, 1) == sp.from("IXI"));
        CHECK(sp.from("IZI").cx(0, 1) == sp.from("ZZI"));
        CHECK(sp.from("ZII").cx(0, 1) == sp.from("ZII"));
        CHECK(sp.from("XII").cz(0, 1) == sp.from("XZI"));
        CHECK(sp.from("XII").cy(0, 1) == sp.from("XYI"));
        CHECK(sp.from("XZI").swap(0, 1) == sp.from("ZXI"));

        // Pauli gates do nothing without phases.
        CHECK(sp.from("XYZ").x(0) == sp.from("XYZ"));
        CHECK(sp.from("XYZ").y(1) == sp.from("XYZ"));
        CHECK(sp.from("XYZ").z(2) == sp.from("XYZ"));
    }

    SECTION("bijection and involution");
    {
        PauliSpace sp(n);
        const PauliSet b = sp.from_list({"XIZ", "IYI", "ZZZ", "IIX"});
        for (const GateSpec &g : gates) {
            CHECK_AT(g.apply(b).size() == b.size(), (std::string(g.name) + " preserves cardinality").c_str());
            CHECK_AT(g.apply(g.apply(b)) == b, (std::string(g.name) + " is an involution").c_str());
        }
    }

    SECTION("simultaneity");
    {
        PauliSpace sp(n);
        // cx moves x forward and z backward at once. Doing the two assignments
        // one after the other would send X_c to X_c X_t and then pick the new
        // x_t up again; the check that catches it is simply that cx undoes
        // itself, together with the table above.
        CHECK(sp.from("XIZ").cx(0, 1).cx(0, 1) == sp.from("XIZ"));
        CHECK(sp.from("XXI").cx(0, 1) == sp.from("XII"));   // not "XXI" or "XYI"

        // swap is four renamings that all read the original.
        CHECK(sp.from("XZI").swap(0, 1).swap(0, 1) == sp.from("XZI"));
    }

    SECTION("commutation is preserved");
    {
        PauliSpace sp(n);
        // A symplectic map sends the set commuting with P to the set commuting
        // with the image of P, which is why a stabilizer group stays one.
        for (const GateSpec &g : gates) {
            const std::string p       = "XYZ";
            const Strings     image_p = spbdd_ref::conjugate({p}, g.qubits, g.images, n);
            CHECK_AT(g.apply(sp.commuting_with(p)) == sp.commuting_with(*image_p.begin()),
                     (std::string(g.name) + " preserves commutation").c_str());
        }
    }

    SECTION("a circuit on the Steane code");
    {
        PauliSpace                     s7(7);
        const std::vector<std::string> gen = {"IIIXXXX", "IXXIIXX", "XIXIXIX",
                                              "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"};
        const PauliSet S = s7.generated_by(gen);
        const PauliSet N = s7.commuting_with_all(gen);

        // Transversal Hadamard maps the Steane code to itself: it swaps the X
        // and Z stabiliser families, which for this code is a relabelling.
        PauliSet SH = S;
        for (int q = 0; q < 7; ++q) SH = SH.h(q);
        CHECK(SH == S);

        // A qubit permutation preserves weight too, so it also keeps the
        // distance -- it is only a relabelling of the physical qubits.
        PauliSet Sp = S, Np = N;
        Sp = Sp.swap(0, 3).swap(1, 5);
        Np = Np.swap(0, 3).swap(1, 5);
        CHECK((Np - Sp).min_weight() == 3);

        // An entangling circuit is still a symplectic bijection, so everything
        // that depends only on the group structure survives it...
        PauliSet Sc = S, Nc = N;
        const std::vector<std::pair<int, int>> circuit = {{0, 1}, {2, 3}, {4, 5}, {1, 4}, {3, 6}};
        for (const auto &e : circuit) { Sc = Sc.cx(e.first, e.second); Nc = Nc.cx(e.first, e.second); }
        for (int q : {0, 2, 5})       { Sc = Sc.s(q);                  Nc = Nc.s(q); }

        CHECK(Sc.size() == 64.0);
        CHECK(Nc.size() == 256.0);
        CHECK(Sc.subset_of(Nc));
        CHECK(Sc != S);

        // ...but weight is not one of those things. CX spreads support across
        // qubits, so the image is a different code with a different distance --
        // this one falls to 1. Only weight-preserving gates (single-qubit
        // Cliffords and qubit permutations) leave the distance alone.
        CHECK((Nc - Sc).min_weight() == 1);
    }

    SECTION("argument checking");
    {
        PauliSpace sp(n);
        CHECK_THROWS(sp.all().h(9), std::out_of_range);
        CHECK_THROWS(sp.all().x(-1), std::out_of_range);
        CHECK_THROWS(sp.all().cx(0, 0), std::invalid_argument);
        CHECK_THROWS(sp.all().cz(1, 1), std::invalid_argument);
        CHECK_THROWS(sp.all().swap(0, 7), std::out_of_range);
    }

    SECTION("reference counting");
    {
        PauliSpace s4(4);
        {
            PauliSet b = s4.weight_at_most(2);
            for (int q = 0; q < 4; ++q) b = b.h(q).s(q);
            b = b.cx(0, 1).cy(1, 2).cz(2, 3).swap(0, 3);
            (void)b.size();
        }
        CHECK(s4.manager().check_zero_ref() == 0);
    }

    return REPORT("gate_test");
}
