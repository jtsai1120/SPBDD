// ===========================================================================
//  Manager and Bdd: the CUDD wrapper.
//
//  Where smoke_test checks that CUDD itself behaves, this checks that the
//  wrapper around it does -- in particular the two contracts stated in
//  bdd.hpp, reference counting and canonicity.
// ===========================================================================

#include "check.hpp"

#include <spbdd/spbdd.hpp>

#include <stdexcept>
#include <utility>

using namespace spbdd;

int main()
{
    Manager mgr(4);

    {
        Bdd x0 = mgr.literal(0), x1 = mgr.literal(1);
        Bdd f = x0 & x1;

        SECTION("variables");
        CHECK(mgr.var_count() == 4);
        mgr.ensure_var_count(6);
        CHECK(mgr.var_count() == 6);
        CHECK_THROWS(mgr.literal(-1), std::out_of_range);
        CHECK_THROWS(mgr.var_to_level(99), std::out_of_range);

        SECTION("boolean algebra");
        // Two ways of building the same function must land on one node; every
        // set comparison in the library above rests on this.
        CHECK(f == mgr.literal(0).ite(mgr.literal(1), mgr.constant(false)));
        CHECK(f.node_count() == 3);
        CHECK(f.sat_count(4) == 4.0);
        CHECK((f | !f).is_true());
        CHECK((f & !f).is_false());
        CHECK(f.diff(x1).is_false());
        CHECK((x0 ^ x1) == ((x0 | x1) & !(x0 & x1)));

        SECTION("handle lifetime");
        Bdd g = f;
        g     = g;                  // self-assignment must not drop the node
        Bdd h = std::move(g);
        CHECK(h == f);
        CHECK(!Bdd().valid());
        CHECK_THROWS(Bdd() & f, std::logic_error);

        SECTION("structure");
        CHECK(f.top_var() == 0);
        CHECK(f.high() == x1);
        CHECK(f.low().is_false());
        CHECK_THROWS(mgr.constant(true).top_var(), std::logic_error);

        // A complemented parent has to hand complemented children back, which
        // is what the whole enumeration layer depends on.
        Bdd nf = !f;
        CHECK(nf.high() == !x1);
        CHECK(nf.low().is_true());

        SECTION("quantification");
        CHECK(f.exists({1}) == x0);
        CHECK(f.forall({1}).is_false());
        CHECK(f.exists({}) == f);
        CHECK_THROWS(f.exists({99}), std::out_of_range);

        SECTION("simultaneous substitution");
        // {x0 := x0^x1, x1 := x0} applied to (x0 & x1). Done simultaneously
        // the answer is (x0^x1) & x0; done one variable at a time it collapses
        // to false, so this distinguishes the two.
        Bdd sub = f.compose({{0, x0 ^ x1}, {1, x0}});
        CHECK(sub == (x0 & !x1));
        CHECK(!sub.is_false());
        CHECK(f.compose({}) == f);

        // permute renames variables rather than substituting functions
        CHECK(f.permute({{0, 2}, {1, 3}}) == (mgr.literal(2) & mgr.literal(3)));

        SECTION("reordering");
        // Sifting permutes levels, never variable numbers, and leaves
        // referenced nodes valid.
        const double before = f.sat_count(4);
        const int    level0 = mgr.var_to_level(0);
        mgr.reorder_now();
        CHECK(f.sat_count(4) == before);
        CHECK(f == (mgr.literal(0) & mgr.literal(1)));
        CHECK(mgr.level_to_var(mgr.var_to_level(0)) == 0);
        (void)level0;

        SECTION("managers are separate");
        Manager other(4);
        CHECK_THROWS(f & other.literal(0), std::logic_error);
        CHECK_THROWS(f.compose({{0, other.literal(0)}}), std::logic_error);
    }

    SECTION("reference counting");
    // Every Bdd above has now been destroyed, so CUDD should hold nothing.
    CHECK(mgr.check_zero_ref() == 0);

    return REPORT("bdd_test");
}
