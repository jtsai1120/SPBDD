// ===========================================================================
//  smoke_test -- proves the build reaches CUDD, and that the specific CUDD
//  behaviours this library is going to rely on actually hold.
//
//  It talks to CUDD directly and deliberately keeps doing so once the wrapper
//  in spbdd/bdd.hpp exists: when something breaks, this test is the one that
//  says whether the problem is ours or the package's.
// ===========================================================================

// cudd.h uses size_t and FILE but includes neither header itself, so these two
// lines are load-bearing and must come first.
#include <cstddef>
#include <cstdio>

#include <cudd.h>

#include "check.hpp"

int main()
{
    DdManager *m = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    if (!m) {
        std::printf("Cudd_Init failed\n");
        return 1;
    }

    // --- variables ---------------------------------------------------------
    // Four boolean variables: two qubits' worth of x0 z0 x1 z1. CUDD creates a
    // variable the first time it is asked for, so there is no "declare" step.
    std::printf("variables\n");
    const int NV = 4;
    DdNode   *v[NV];
    for (int i = 0; i < NV; ++i) {
        v[i] = Cudd_bddIthVar(m, i);
        Cudd_Ref(v[i]);
    }
    CHECK(Cudd_ReadSize(m) == NV);
    CHECK(Cudd_NodeReadIndex(v[2]) == 2);

    // --- boolean algebra and canonicity ------------------------------------
    // Two different ways of building the same function must land on the same
    // node. Everything downstream (set equality as a pointer compare) rests on
    // this, so it is worth asserting once rather than assuming.
    std::printf("boolean algebra\n");
    DdNode *f = Cudd_bddAnd(m, v[0], v[1]);
    Cudd_Ref(f);
    DdNode *g = Cudd_bddIte(m, v[0], v[1], Cudd_ReadLogicZero(m));
    Cudd_Ref(g);
    CHECK(f == g);
    CHECK(Cudd_DagSize(f) == 3);
    CHECK(Cudd_CountMinterm(m, f, NV) == 4.0);   // x0 x1 fixed, x2 x3 free

    // --- complement edges --------------------------------------------------
    // Negation is a bit flipped in the pointer, not a new node. Any traversal
    // we write has to resolve this or it will read the wrong cofactors.
    std::printf("complement edges\n");
    DdNode *nf = Cudd_Not(f);
    CHECK(Cudd_Regular(nf) == Cudd_Regular(f));
    CHECK(Cudd_IsComplement(nf) != Cudd_IsComplement(f));

    // --- quantification ----------------------------------------------------
    // A cube is itself a BDD (the conjunction of the variables to remove), and
    // it needs its own reference like anything else.
    std::printf("quantification\n");
    DdNode *cube = v[1];
    Cudd_Ref(cube);
    DdNode *ex = Cudd_bddExistAbstract(m, f, cube);
    Cudd_Ref(ex);
    CHECK(ex == v[0]);                            // exists x1 . (x0 & x1) == x0

    // --- reordering --------------------------------------------------------
    // Sifting rewrites levels, never variable indices, and leaves referenced
    // nodes valid. Code that indexes by variable number therefore survives it;
    // code that assumes level == index does not, which is what this catches.
    std::printf("reordering\n");
    const double minterms_before = Cudd_CountMinterm(m, f, NV);
    Cudd_ReduceHeap(m, CUDD_REORDER_SIFT, 0);
    CHECK(Cudd_CountMinterm(m, f, NV) == minterms_before);
    CHECK(Cudd_NodeReadIndex(v[0]) == 0);
    DdNode *f2 = Cudd_bddAnd(m, v[0], v[1]);
    Cudd_Ref(f2);
    CHECK(f2 == f);

    // --- reference counting ------------------------------------------------
    // Hand every reference back, then let CUDD audit us. A non-zero count here
    // is a leak, and it is the check the Bdd RAII wrapper will have to keep
    // passing once it exists.
    std::printf("reference counting\n");
    Cudd_RecursiveDeref(m, f2);
    Cudd_RecursiveDeref(m, ex);
    Cudd_RecursiveDeref(m, cube);
    Cudd_RecursiveDeref(m, g);
    Cudd_RecursiveDeref(m, f);
    for (int i = 0; i < NV; ++i) Cudd_RecursiveDeref(m, v[i]);
    CHECK(Cudd_CheckZeroRef(m) == 0);

    Cudd_Quit(m);

    return REPORT("smoke_test");
}
