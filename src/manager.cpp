// ===========================================================================
//  Implementation of spbdd::Manager.
// ===========================================================================

#include "spbdd/manager.hpp"

// cudd.h uses size_t and FILE but includes neither header itself, so these two
// lines are load-bearing and must come first.
#include <cstddef>
#include <cstdio>

#include <cudd.h>

#include <stdexcept>
#include <string>

namespace spbdd {
namespace {

// CUDD signals failure by returning null (typically out of memory, or a
// reordering that was interrupted). Turning that into an exception here keeps
// the null out of every call site.
DdNode *checked(DdNode *node, const char *what)
{
    if (!node) throw std::runtime_error(std::string("CUDD failed in ") + what);
    return node;
}

} // namespace

Manager::Manager(int n_vars, const ManagerConfig &cfg)
{
    if (n_vars < 0) throw std::invalid_argument("Manager: negative variable count");

    dd_ = Cudd_Init(static_cast<unsigned>(n_vars),
                    0,                                                       // no ZDD variables
                    cfg.unique_slots ? cfg.unique_slots : CUDD_UNIQUE_SLOTS,
                    cfg.cache_size ? cfg.cache_size : CUDD_CACHE_SLOTS,
                    cfg.max_memory);
    if (!dd_) throw std::runtime_error("Cudd_Init failed");

    if (cfg.dynamic_reordering) Cudd_AutodynEnable(dd_, CUDD_REORDER_SIFT);
}

Manager::~Manager()
{
    if (dd_) Cudd_Quit(dd_);
}

int Manager::var_count() const
{
    return Cudd_ReadSize(dd_);
}

void Manager::ensure_var_count(int n_vars)
{
    if (n_vars < 0) throw std::invalid_argument("ensure_var_count: negative variable count");
    // Asking for a variable creates it and everything below it. Shrinking is
    // not offered because CUDD cannot do it safely.
    if (n_vars > Cudd_ReadSize(dd_))
        checked(Cudd_bddIthVar(dd_, n_vars - 1), "Cudd_bddIthVar");
}

Bdd Manager::constant(bool value)
{
    DdNode *n = value ? Cudd_ReadOne(dd_) : Cudd_ReadLogicZero(dd_);
    Cudd_Ref(n);
    return Bdd::adopt(this, n);
}

Bdd Manager::literal(int var, bool positive)
{
    if (var < 0) throw std::out_of_range("literal: negative variable number");

    DdNode *n = checked(Cudd_bddIthVar(dd_, var), "Cudd_bddIthVar");
    if (!positive) n = Cudd_Not(n);   // negation is a mark on the edge, not a node
    Cudd_Ref(n);
    return Bdd::adopt(this, n);
}

void Manager::set_dynamic_reordering(bool enabled)
{
    if (enabled) Cudd_AutodynEnable(dd_, CUDD_REORDER_SIFT);
    else         Cudd_AutodynDisable(dd_);
}

void Manager::reorder_now()
{
    Cudd_ReduceHeap(dd_, CUDD_REORDER_SIFT, 0);
}

int Manager::var_to_level(int var) const
{
    if (var < 0 || var >= Cudd_ReadSize(dd_))
        throw std::out_of_range("var_to_level: variable outside the declared space");
    return Cudd_ReadPerm(dd_, var);
}

int Manager::level_to_var(int level) const
{
    if (level < 0 || level >= Cudd_ReadSize(dd_))
        throw std::out_of_range("level_to_var: level outside the declared space");
    return Cudd_ReadInvPerm(dd_, level);
}

std::size_t Manager::live_nodes() const
{
    return static_cast<std::size_t>(Cudd_ReadNodeCount(dd_));
}

std::size_t Manager::peak_nodes() const
{
    return static_cast<std::size_t>(Cudd_ReadPeakNodeCount(dd_));
}

std::size_t Manager::memory_in_use() const
{
    return Cudd_ReadMemoryInUse(dd_);
}

int Manager::check_zero_ref() const
{
    return Cudd_CheckZeroRef(dd_);
}

} // namespace spbdd
