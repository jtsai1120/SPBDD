// ===========================================================================
//  spbdd by example: the Steane [[7,1,3]] code, end to end.
//
//  Builds against the public headers and libspbdd.a alone:
//
//      g++ -std=c++17 -Iinclude examples/demo.cpp build/libspbdd.a -lm -o demo
//
//  Nothing here mentions CUDD, and nothing here enumerates a set except where
//  the point is to print its elements.
// ===========================================================================

#include <spbdd/spbdd.hpp>

#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace spbdd;

namespace {

void rule(const char *title)
{
    std::printf("\n%s\n", title);
    for (const char *c = title; *c; ++c) std::printf("-");
    std::printf("\n");
}

// Cardinality and diagram size side by side: the premise of the whole library
// is that these two numbers are unrelated.
void describe(const char *name, const PauliSet &s)
{
    std::printf("  %-28s %20.0f elements  %5zu nodes\n", name, s.size(), s.node_count());
}

std::string bits(const std::vector<bool> &v)
{
    std::string s;
    for (bool b : v) s.push_back(b ? '1' : '0');
    return s;
}

} // namespace

static void steane_walkthrough()
{
    std::printf("spbdd -- the Steane [[7,1,3]] code\n");

    const std::vector<std::string> generators = {"IIIXXXX", "IXXIIXX", "XIXIXIX",
                                                 "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"};

    PauliSpace           sp(7);
    const StabilizerCode code(sp, generators);

    // -----------------------------------------------------------------------
    rule("code data (computed once)");

    std::printf("  %d qubits, %d stabilizers, %d logical qubit(s), distance %d\n\n",
                code.n_qubits(), code.n_stabilizers(), code.n_logical(), code.distance());

    std::printf("  %-14s %-14s\n", "stabilizer", "destabilizer");
    for (int i = 0; i < code.n_stabilizers(); ++i)
        std::printf("  %-14s %-14s\n", code.stabilizers()[static_cast<std::size_t>(i)].c_str(),
                    code.destabilizers()[static_cast<std::size_t>(i)].c_str());
    std::printf("\n  logical X %s      logical Z %s\n\n", code.logical_x()[0].c_str(),
                code.logical_z()[0].c_str());

    describe("S (the stabilizer group)", code.group());
    describe("N(S) (its normaliser)", code.normalizer());
    describe("N(S) - S (logical ops)", code.normalizer() - code.group());

    // The distance is the lightest logical operator, found by bisecting on the
    // weight bound: three probes here, and no element is ever looked at.
    std::printf("\n  distance = (N(S) - S).min_weight() = %d\n",
                (code.normalizer() - code.group()).min_weight());

    // -----------------------------------------------------------------------
    rule("cardinality and diagram size are unrelated");

    describe("all Paulis on 7 qubits", sp.all());
    describe("weight <= 1", sp.weight_at_most(1));
    describe("weight <= 3", sp.weight_at_most(3));
    describe("X on qubit 0, rest free", sp.matching("X******"));

    // -----------------------------------------------------------------------
    rule("set algebra");

    const PauliSet single = sp.weight_at_most(1);
    describe("single-qubit errors", single);

    std::printf("  the ones on qubit 0: ");
    for (const std::string &s : single.to_strings())
        if (s[0] != 'I') std::printf("%s ", s.c_str());
    std::printf("\n");

    // Distance 3 means exactly this: no two single-qubit errors multiply into a
    // logical operator, so a decoder can always tell them apart.
    std::printf("  are any two of them indistinguishable? %s\n",
                code.has_inequivalent_pair(single) ? "yes" : "no");

    // -----------------------------------------------------------------------
    rule("gates act by conjugation");

    PauliSet transversal_h = code.group();
    for (int q = 0; q < 7; ++q) transversal_h = transversal_h.h(q);
    std::printf("  transversal H maps S to %s\n",
                transversal_h == code.group() ? "itself" : "a different group");

    // Conjugating the generators one at a time gives the generators of another
    // code. The circuit is a bijection, so the group keeps its size -- but CX
    // moves support between qubits, and weight is not preserved.
    const std::vector<std::pair<int, int>> circuit = {{0, 1}, {2, 3}, {4, 5}, {1, 4}};
    std::vector<std::string>               moved_generators;
    for (const std::string &g : generators) {
        PauliSet one = sp.from(g);
        for (const auto &e : circuit) one = one.cx(e.first, e.second);
        moved_generators.push_back(*one.any_element());
    }
    const StabilizerCode moved(sp, moved_generators);
    std::printf("  a CX circuit gives another code: |S| = %.0f still, distance %d -> %d\n",
                moved.group().size(), code.distance(), moved.distance());

    // -----------------------------------------------------------------------
    rule("circuit level");

    const PauliSet clean = sp.from_list({"XIIIIII", "IIZIIII"});
    describe("two known errors", clean);

    // A fault at a two-qubit gate location multiplies by all 16 Paulis on those
    // qubits. That is one existential quantification rather than 16
    // multiplications -- note the element count rising while the diagram
    // shrinks, because those variables stop being asked about at all.
    const PauliSet faulted = clean.forget({2, 3});
    describe("after a fault on {2,3}", faulted);

    // Measuring splits the set rather than changing it. Measuring a stabilizer
    // generator is one bit of syndrome extraction.
    const MeasurementSplit m = faulted.measure(generators[0]);
    describe("  outcome as expected", m.unflipped);
    describe("  outcome flipped", m.flipped);

    // Reset is quantify-then-pin, in that order: a real projection, so the set
    // shrinks back.
    describe("after reset of {2,3}", faulted.reset({2, 3}));

    // -----------------------------------------------------------------------
    rule("is a decoder safe on this error set?");

    for (int w : {1, 2}) {
        const PauliSet errors = sp.weight_at_most(w);
        const auto     found  = code.find_inequivalent_pair(errors);

        std::printf("  weight <= %d  (%.0f errors, %zu nodes):  ", w, errors.size(),
                    errors.node_count());
        if (!found) {
            std::printf("safe\n");
            continue;
        }

        const std::string product = pauli_string_to_text(
            pauli_mul(parse_pauli_string(found->first), parse_pauli_string(found->second)));
        std::printf("UNSAFE\n");
        std::printf("      %s and %s share the syndrome %s,\n", found->first.c_str(),
                    found->second.c_str(), bits(found->syndrome).c_str());
        std::printf("      but their logical signatures differ: %s vs %s.\n",
                    bits(found->signature_first).c_str(), bits(found->signature_second).c_str());
        std::printf("      Their product %s is a logical operator, so no decoder can\n",
                    product.c_str());
        std::printf("      tell them apart -- yet correcting one leaves a logical error.\n");
    }

    std::printf("\n");
}

// A second run, with its own space: this backend allows only one at a time.
static void scale()
{
    rule("the extreme case");

    PauliSpace big(30);
    describe("all Paulis on 30 qubits", big.all());
    describe("...of weight <= 2", big.weight_at_most(2));
    std::printf("\n");
}

int main()
{
    steane_walkthrough();
    scale();
    return 0;
}
