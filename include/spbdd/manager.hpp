#pragma once

// ===========================================================================
//  spbdd::Manager -- the BuDDy instance, and the factory for the Bdds living
//  in it.
//
//  BuDDy's own header is deliberately absent here, for the reason given in
//  bdd.hpp, which also states the reference-counting and canonicity contracts
//  these two classes maintain together.
//
//  *** BuDDy keeps all of its state in process-wide globals. ***
//
//  There is no manager object to hold: bdd_init() takes no handle and there is
//  one variable table for the whole program. So at most one Manager -- and
//  therefore at most one PauliSpace -- may be alive at a time, and building a
//  second one throws rather than quietly resetting the first one's table out
//  from under it. Code that needs two universes has to finish with one before
//  starting the other.
// ===========================================================================

#include "spbdd/bdd.hpp"

#include <cstddef>

namespace spbdd {

// Start-up hints. Every field is a suggestion: BuDDy grows its table on
// demand, so these only trade start-up memory against the number of rehashes
// a long run pays for.
struct ManagerConfig {
    // 0 means "a default that suits a small problem".
    unsigned    initial_nodes = 0;   // size of the node table at bdd_init
    unsigned    cache_size    = 0;   // operator cache size

    // Ignored by this backend: BuDDy caps the node count rather than the byte
    // count, and has no equivalent knob. Kept so the field survives a switch
    // of backend.
    std::size_t max_memory = 0;

    // Sifting. Off by default: it is a large win on some problems, but it
    // makes node counts depend on when a reorder happened, which makes
    // measurements harder to reproduce, and on the workloads measured so far
    // it costs more time than it saves. Turn it on deliberately.
    bool dynamic_reordering = false;
};

class Manager {
public:
    // Throws std::runtime_error if another Manager is already alive.
    explicit Manager(int n_vars = 0, const ManagerConfig &cfg = {});
    ~Manager();

    Manager(const Manager &)            = delete;
    Manager &operator=(const Manager &) = delete;

    // True while some Manager exists, which is the condition that makes
    // constructing another one throw.
    static bool alive() noexcept;

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
    // code that confuses the two silently stops working. BuDDy offers only
    // level-permuting methods, so unlike CUDD there is no method here that can
    // rewrite what a variable means.
    void set_dynamic_reordering(bool enabled);
    void reorder_now();
    int  var_to_level(int var) const;
    int  level_to_var(int level) const;

    // --- statistics --------------------------------------------------------
    std::size_t live_nodes() const;
    std::size_t peak_nodes() const;
    std::size_t memory_in_use() const;

    // Nodes still referenced beyond the ones BuDDy keeps permanently for the
    // variable projections. Zero once every Bdd is gone, which makes this the
    // acceptance test for the Bdd class. Unlike CUDD's Cudd_CheckZeroRef this
    // is a live-node count against a baseline rather than a reference-count
    // audit, so it catches leaks but cannot attribute them.
    int check_zero_ref() const;

private:
    std::size_t baseline_nodes_ = 0;   // live nodes owned by the variable table
    std::size_t peak_nodes_     = 0;
};

} // namespace spbdd
