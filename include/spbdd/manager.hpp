#pragma once

// ===========================================================================
//  spbdd::Manager -- one CUDD instance, and the factory for the Bdds living
//  in it.
//
//  CUDD's own header is deliberately absent here; DdManager is forward-
//  declared for the same reason as DdNode in bdd.hpp, which also states the
//  reference-counting and canonicity contracts these two classes maintain
//  together.
//
//  Unlike some BDD packages CUDD keeps no process-wide state, so several
//  Managers may be alive at once. Functions from different Managers share no
//  unique table and may never be mixed; every operation in Bdd checks it.
// ===========================================================================

#include "spbdd/bdd.hpp"

#include <cstddef>

struct DdManager;   // CUDD's manager type, forward-declared on purpose

namespace spbdd {

// Start-up hints. Every field is a suggestion: CUDD grows its tables on
// demand, so these only trade start-up memory against the number of rehashes
// a long run pays for.
struct ManagerConfig {
    // 0 means "whatever CUDD's own default is".
    unsigned    unique_slots = 0;   // initial size of each variable's subtable
    unsigned    cache_size   = 0;   // initial operator-cache size
    std::size_t max_memory   = 0;   // soft cap in bytes

    // Sifting. Off by default: it is a large win on some problems, but it
    // makes node counts depend on when a reorder happened, which makes
    // measurements harder to reproduce. Turn it on deliberately.
    bool dynamic_reordering = false;
};

class Manager {
public:
    explicit Manager(int n_vars = 0, const ManagerConfig &cfg = {});
    ~Manager();

    Manager(const Manager &)            = delete;
    Manager &operator=(const Manager &) = delete;

    // --- variables ---------------------------------------------------------
    // Variables are plain non-negative integers. The space can only grow: no
    // package can take a variable back while live functions may mention it.
    int  var_count() const;
    void ensure_var_count(int n_vars);

    // --- construction ------------------------------------------------------
    Bdd constant(bool value);
    Bdd literal(int var, bool positive = true);

    // --- variable order ----------------------------------------------------
    // Reordering permutes *levels*. Variable numbers never change, so code
    // written against variable numbers keeps working across a reorder -- and
    // code that confuses the two silently stops working.
    void set_dynamic_reordering(bool enabled);
    void reorder_now();
    int  var_to_level(int var) const;
    int  level_to_var(int level) const;

    // --- statistics --------------------------------------------------------
    std::size_t live_nodes() const;
    std::size_t peak_nodes() const;
    std::size_t memory_in_use() const;

    // Nodes still referenced. Zero once every Bdd is gone; anything else is a
    // leak in this library rather than in the caller's code, which makes this
    // the acceptance test for the Bdd class.
    int check_zero_ref() const;

    DdManager *raw() const noexcept { return dd_; }   // escape hatch

private:
    DdManager *dd_ = nullptr;
};

} // namespace spbdd
