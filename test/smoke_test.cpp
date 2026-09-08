// ===========================================================================
//  smoke_test -- proves the build reaches BuDDy, and that the specific BuDDy
//  behaviours this library relies on actually hold.
//
//  It talks to BuDDy directly and deliberately keeps doing so alongside the
//  wrapper in spbdd/bdd.hpp: when something breaks, this test is the one that
//  says whether the problem is ours or the package's.
// ===========================================================================

#include <bdd.h>

// BuDDy's header, compiled as C++, redirects part of its C API to overloads
// returning objects of its own bdd class. This test wants the plain C API, the
// same one the library is written against.
#undef bddtrue
#undef bddfalse
#undef bdd_init
#undef bdd_ithvar
#undef bdd_nithvar
#undef bdd_makeset
#undef bdd_ibuildcube
#undef bdd_anodecount
extern "C" const BDD bddtrue;
extern "C" const BDD bddfalse;

#include "check.hpp"

namespace {

BDD claim(BDD r)
{
    bdd_addref(r);
    return r;
}

void quiet_gbc(int, bddGbcStat *) {}

} // namespace

int main()
{
    if (bdd_init(100000, 10000) < 0) {
        std::printf("bdd_init failed\n");
        return 1;
    }
    bdd_gbc_hook(quiet_gbc);

    // --- variables ---------------------------------------------------------
    // Unlike CUDD, BuDDy needs its variables declared up front; asking for one
    // beyond the declared count is an error rather than a request to grow.
    SECTION("variables");
    const int NV = 4;
    bdd_setvarnum(NV);
    CHECK(bdd_varnum() == NV);
    CHECK(bdd_var(bdd_ithvar(2)) == 2);

    // --- boolean algebra and canonicity ------------------------------------
    // Two different ways of building the same function must land on the same
    // node. Everything downstream (set equality as an integer compare) rests
    // on this, so it is worth asserting once rather than assuming.
    SECTION("boolean algebra");
    BDD x0 = claim(bdd_ithvar(0));
    BDD x1 = claim(bdd_ithvar(1));
    BDD f  = claim(bdd_and(x0, x1));
    BDD g  = claim(bdd_ite(x0, x1, bddfalse));
    CHECK(f == g);
    CHECK(bdd_nodecount(f) == 2);            // internal nodes only; no terminal
    CHECK(bdd_satcount(f) == 4.0);           // over all four declared variables

    // --- no complement edges -----------------------------------------------
    // BuDDy stores negation as a distinct node rather than a mark on an edge,
    // so a traversal reads the plain Shannon decomposition with nothing to
    // resolve -- which is why Bdd::low/high have no mark handling.
    SECTION("no complement edges");
    BDD nf = claim(bdd_not(f));
    CHECK(nf != f);
    CHECK(bdd_high(f) == x1);
    CHECK(bdd_low(f) == bddfalse);

    // --- quantification ----------------------------------------------------
    // A varset is itself a BDD, and it needs its own reference like anything
    // else.
    SECTION("quantification");
    BDD cube = claim(bdd_ithvar(1));
    BDD ex   = claim(bdd_exist(f, cube));
    CHECK(ex == x0);                          // exists x1 . (x0 & x1) == x0

    // --- simultaneous substitution -----------------------------------------
    // A fresh pair is the identity on every variable, and veccompose reads
    // every right-hand side in the original function. Applying
    // {x0 := x0^x1, x1 := x0} to (x0 & x1) gives x0 & !x1; one assignment at a
    // time would collapse it to false, so this distinguishes the two.
    SECTION("simultaneous substitution");
    bddPair *pair  = bdd_newpair();
    BDD      xor01 = claim(bdd_xor(x0, x1));
    bdd_setbddpair(pair, 0, xor01);
    bdd_setbddpair(pair, 1, x0);
    BDD sub  = claim(bdd_veccompose(f, pair));
    BDD want = claim(bdd_and(x0, bdd_nithvar(1)));
    CHECK(sub == want);
    CHECK(sub != bddfalse);

    // A renaming that reverses the order of two variables has to come out as a
    // transposition rather than collapsing them, which is what swap() needs.
    bddPair *swap = bdd_newpair();
    bdd_setpair(swap, 0, 1);
    bdd_setpair(swap, 1, 0);
    BDD swapped  = claim(bdd_replace(f, swap));
    CHECK(swapped == f);                       // x0 & x1 is symmetric
    BDD asym     = claim(bdd_and(x0, bdd_nithvar(1)));
    BDD asym_sw  = claim(bdd_replace(asym, swap));
    BDD asym_want = claim(bdd_and(x1, bdd_nithvar(0)));
    CHECK(asym_sw == asym_want);

    // --- reordering --------------------------------------------------------
    // Sifting rewrites levels, never variable numbers, and leaves referenced
    // nodes valid. BuDDy needs variable blocks declared before it will move
    // anything at all, which is why the Manager declares them at start-up.
    SECTION("reordering");
    const double minterms_before = bdd_satcount(f);
    bdd_varblockall();
    bdd_reorder(BDD_REORDER_SIFT);
    CHECK(bdd_satcount(f) == minterms_before);
    CHECK(bdd_var(bdd_ithvar(0)) == 0);
    BDD f2 = claim(bdd_and(bdd_ithvar(0), bdd_ithvar(1)));
    CHECK(f2 == f);

    // --- reference counting ------------------------------------------------
    // Hand every reference back, then let a collection show the nodes go away.
    // BuDDy has no per-node audit like Cudd_CheckZeroRef, so this is a live
    // count rather than a reference check -- the same substitute Manager uses.
    SECTION("reference counting");
    bdd_gbc();
    const int live_while_held = bdd_getnodenum();

    bdd_delref(asym_want); bdd_delref(asym_sw); bdd_delref(asym);
    bdd_delref(swapped);   bdd_delref(want);    bdd_delref(sub);
    bdd_delref(xor01);     bdd_delref(ex);      bdd_delref(cube);
    bdd_delref(nf);        bdd_delref(f2);      bdd_delref(g);
    bdd_delref(f);         bdd_delref(x1);      bdd_delref(x0);
    bdd_freepair(pair);    bdd_freepair(swap);

    bdd_gbc();
    CHECK(bdd_getnodenum() < live_while_held);

    bdd_done();
    return REPORT("smoke_test");
}
