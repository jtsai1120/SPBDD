#pragma once

// ===========================================================================
//  A test harness small enough to read in one sitting.
//
//  Each .cpp under test/ is built into its own binary and run by `make check`,
//  which looks only at the exit status. So the job here is to count failures,
//  print enough to find one, and return non-zero.
//
//      SECTION("weight");
//      CHECK(sp.weight_at_most(1).size() == 13.0);
//      CHECK_AT(a == b, "weight_exactly(3)");        // inside a loop
//      CHECK_THROWS(s.to_strings(3), std::length_error);
//      return REPORT("pauliset_test");
// ===========================================================================

#include <cstdio>

namespace spbdd_test {

inline int failures = 0;

inline void section(const char *name)
{
    std::printf("%s\n", name);
}

inline void pass(const char *what)
{
    std::printf("  ok    %s\n", what);
}

inline void fail(const char *what, const char *file, int line)
{
    std::printf("  FAIL  %s   (%s:%d)\n", what, file, line);
    ++failures;
}

inline int report(const char *name)
{
    if (failures) {
        std::printf("\n%s: %d check(s) FAILED\n", name, failures);
        return 1;
    }
    std::printf("\n%s: all checks passed\n", name);
    return 0;
}

} // namespace spbdd_test

#define SECTION(name) spbdd_test::section(name)
#define REPORT(name)  spbdd_test::report(name)

// The condition doubles as its own description, which is why the checks in
// these tests are written to read as statements about the library.
#define CHECK(cond)                                                            \
    do {                                                                       \
        if (cond) spbdd_test::pass(#cond);                                     \
        else      spbdd_test::fail(#cond, __FILE__, __LINE__);                 \
    } while (0)

// For checks inside a loop, where the expression is the same every time and
// the interesting part is which iteration failed.
#define CHECK_AT(cond, label)                                                  \
    do {                                                                       \
        if (cond) spbdd_test::pass(label);                                      \
        else      spbdd_test::fail(label, __FILE__, __LINE__);                  \
    } while (0)

#define CHECK_THROWS(expr, exception_type)                                     \
    do {                                                                       \
        bool spbdd_threw_ = false;                                             \
        try { (void)(expr); } catch (const exception_type &) { spbdd_threw_ = true; } \
        if (spbdd_threw_) spbdd_test::pass(#expr " throws " #exception_type);  \
        else spbdd_test::fail(#expr " throws " #exception_type, __FILE__, __LINE__); \
    } while (0)
