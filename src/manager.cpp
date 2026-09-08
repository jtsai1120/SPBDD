// ===========================================================================
//  Implementation of spbdd::Manager against BuDDy.
// ===========================================================================

#include "spbdd/manager.hpp"

#include <bdd.h>

// BuDDy's header, compiled as C++, redirects part of its C API to overloads
// returning objects of its own bdd class -- including the terminals and the
// variable constructors this layer uses. The C entry points are still in the
// library, so take them back: everything below wants plain node indices.
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

#include <stdexcept>
#include <string>

namespace spbdd {
namespace {

// BuDDy's state is process-wide, so this is the whole of the "which manager"
// question: there is at most one.
Manager *g_instance = nullptr;

// Enough for a few hundred thousand nodes without a rehash, which covers the
// codes this library is normally pointed at.
constexpr int kDefaultNodes = 1 << 20;
constexpr int kDefaultCache = 1 << 16;

// BuDDy prints a line on every garbage collection and on every reordering.
// A library has no business writing to the user's terminal.
void silent_gbc(int, bddGbcStat *) {}
void silent_reorder(int) {}

// The live-node count owned by the variable table itself: BuDDy keeps the
// projection function of every variable permanently referenced.
std::size_t table_baseline()
{
    bdd_gbc();
    return static_cast<std::size_t>(bdd_getnodenum());
}

} // namespace

bool Manager::alive() noexcept { return g_instance != nullptr; }

Manager::Manager(int n_vars, const ManagerConfig &cfg)
{
    if (n_vars < 0) throw std::invalid_argument("Manager: negative variable count");
    if (g_instance)
        throw std::runtime_error(
            "Manager: BuDDy keeps its state in process-wide globals, so only one "
            "Manager (and so only one PauliSpace) may be alive at a time");

    const int rc = bdd_init(cfg.initial_nodes ? static_cast<int>(cfg.initial_nodes) : kDefaultNodes,
                            cfg.cache_size ? static_cast<int>(cfg.cache_size) : kDefaultCache);
    if (rc < 0) throw std::runtime_error("bdd_init failed");

    bdd_gbc_hook(silent_gbc);
    bdd_reorder_hook(silent_reorder);

    bdd_setvarnum(n_vars);

    // BuDDy reorders within declared variable blocks and does nothing at all
    // if none exist. One block per variable is the "everything may move" case,
    // which is what the CUDD backend's behaviour corresponds to.
    bdd_varblockall();

    baseline_nodes_ = table_baseline();
    g_instance      = this;
}

Manager::~Manager()
{
    if (g_instance == this) {
        bdd_done();
        g_instance = nullptr;
    }
}

int Manager::var_count() const { return bdd_varnum(); }

void Manager::ensure_var_count(int n_vars)
{
    if (n_vars < 0) throw std::invalid_argument("ensure_var_count: negative variable count");
    if (n_vars <= bdd_varnum()) return;

    bdd_setvarnum(n_vars);
    bdd_varblockall();

    // The new variables bring permanently-live projection nodes with them, so
    // the leak baseline moves too.
    baseline_nodes_ = table_baseline();
}

Bdd Manager::constant(bool value)
{
    const NodeId n = value ? bddtrue : bddfalse;
    bdd_addref(n);
    return Bdd::adopt(this, n);
}

Bdd Manager::literal(int var, bool positive)
{
    if (var < 0) throw std::out_of_range("literal: negative variable number");
    if (var >= bdd_varnum()) ensure_var_count(var + 1);

    const NodeId n = positive ? bdd_ithvar(var) : bdd_nithvar(var);
    bdd_addref(n);
    return Bdd::adopt(this, n);
}

void Manager::set_dynamic_reordering(bool enabled)
{
    if (enabled) bdd_autoreorder(BDD_REORDER_SIFT);
    else         bdd_autoreorder(BDD_REORDER_NONE);
}

void Manager::reorder_now() { bdd_reorder(BDD_REORDER_SIFT); }

int Manager::var_to_level(int var) const
{
    if (var < 0 || var >= bdd_varnum())
        throw std::out_of_range("var_to_level: variable outside the declared space");
    return bdd_var2level(var);
}

int Manager::level_to_var(int level) const
{
    if (level < 0 || level >= bdd_varnum())
        throw std::out_of_range("level_to_var: level outside the declared space");
    return bdd_level2var(level);
}

std::size_t Manager::live_nodes() const
{
    return static_cast<std::size_t>(bdd_getnodenum());
}

std::size_t Manager::peak_nodes() const
{
    // BuDDy reports the allocated table size, which is the high-water mark of
    // what it has had to hold.
    return static_cast<std::size_t>(bdd_getallocnum());
}

std::size_t Manager::memory_in_use() const
{
    // BuDDy has no byte-level accounting; a node is five ints in its table.
    return static_cast<std::size_t>(bdd_getallocnum()) * 5u * sizeof(int);
}

int Manager::check_zero_ref() const
{
    bdd_gbc();
    const std::size_t live = static_cast<std::size_t>(bdd_getnodenum());
    return live > baseline_nodes_ ? static_cast<int>(live - baseline_nodes_) : 0;
}

} // namespace spbdd
