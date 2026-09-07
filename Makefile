# ===========================================================================
#  spbdd -- a BDD-based library for sets of symplectic Pauli operators.
#
#  Layout
#    include/spbdd/   public headers   (installed)
#    src/             implementation   -> build/libspbdd.a
#    test/            one binary per .cpp, run by `make check`
#    examples/        one binary per .cpp, uses the public API only
#    cudd/            CUDD source tree, built in place (not part of this repo)
#
#  Quick start
#    make cudd        fetch + build CUDD into cudd/           (once)
#    make             build build/libspbdd.a
#    make check       build and run every test
#    make config      show what this build is actually using
# ===========================================================================

CXX      ?= g++
AR       ?= ar
CXXSTD   ?= -std=c++17
WARN     ?= -Wall -Wextra
OPT      ?= -O2
CXXFLAGS ?= $(CXXSTD) $(WARN) $(OPT)

BUILD_DIR ?= build
OBJ_DIR   := $(BUILD_DIR)/obj
PREFIX    ?= /usr/local
NPROC     := $(shell nproc 2>/dev/null || echo 1)

# --- CUDD ------------------------------------------------------------------
# Built in place under cudd/ -- no `make install`, nothing lands in /usr/local.
# To use a system-wide CUDD instead:
#     make CUDD_INC=/usr/local/include CUDD_LIB=/usr/local/lib/libcudd.a
CUDD_REPO  ?= https://github.com/ivmai/cudd.git
CUDD_DIR   ?= cudd
CUDD_INC   ?= $(CUDD_DIR)/cudd
CUDD_LIB   ?= $(CUDD_DIR)/cudd/.libs/libcudd.a
CUDD_FLAGS ?= --enable-obj

# --- sources ---------------------------------------------------------------
SRCS := $(wildcard src/*.cpp)
OBJS := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
LIB  := $(BUILD_DIR)/libspbdd.a

TEST_SRCS := $(wildcard test/*.cpp)
TEST_BINS := $(patsubst test/%.cpp,$(BUILD_DIR)/%,$(TEST_SRCS))
EX_SRCS   := $(wildcard examples/*.cpp)
EX_BINS   := $(patsubst examples/%.cpp,$(BUILD_DIR)/%,$(EX_SRCS))

INCLUDES := -Iinclude -I$(CUDD_INC)
LDLIBS   := $(LIB) -lm
DEPFLAGS  = -MMD -MP

ifeq ($(wildcard $(CUDD_LIB)),)
$(warning CUDD not found at $(CUDD_LIB) -- run `make cudd` first)
endif

# ===========================================================================

.PHONY: all lib cudd deps check examples config install uninstall clean distclean help
.DEFAULT_GOAL := all

all: lib

lib: $(LIB)

# libspbdd.a bundles CUDD's objects as well, so a program that uses this
# library links one archive and needs no CUDD on its include or library path.
$(LIB): $(OBJS) $(CUDD_LIB)
	@mkdir -p $(dir $@) $(BUILD_DIR)/cudd_obj
	@rm -f $@ $(BUILD_DIR)/cudd_obj/*.o
	cd $(BUILD_DIR)/cudd_obj && $(AR) x $(abspath $(CUDD_LIB))
	$(AR) rcs $@ $(OBJS) $(BUILD_DIR)/cudd_obj/*.o
	@echo "built $@ ($$($(AR) t $@ | wc -l) objects)"

$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Tests and examples: one binary per source file.
$(BUILD_DIR)/%: test/%.cpp $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Itest $< $(LDLIBS) -o $@

$(BUILD_DIR)/%: examples/%.cpp $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< $(LDLIBS) -o $@

check: $(TEST_BINS)
	@fail=0; \
	for t in $(TEST_BINS); do \
	    printf '%-40s' "$$t"; \
	    if ./$$t >/dev/null 2>&1; then echo "ok"; else echo "FAIL"; fail=1; fi; \
	done; \
	exit $$fail

examples: $(EX_BINS)

# --- CUDD ------------------------------------------------------------------
# Fetch and build CUDD in place. Each step is skipped if already done, so this
# is safe to re-run; `make distclean` throws the build away but keeps the tree.
cudd deps: $(CUDD_LIB)

$(CUDD_LIB):
	@for t in git autoreconf automake libtoolize; do \
	    command -v $$t >/dev/null 2>&1 || { \
	        echo "missing build tool: $$t"; \
	        echo "  sudo apt install -y build-essential git autoconf automake libtool"; \
	        exit 1; }; \
	done
	@test -d $(CUDD_DIR) || git clone --depth 1 $(CUDD_REPO) $(CUDD_DIR)
	@test -f $(CUDD_DIR)/configure || (cd $(CUDD_DIR) && autoreconf -i)
	@test -f $(CUDD_DIR)/config.h  || (cd $(CUDD_DIR) && ./configure $(CUDD_FLAGS))
	$(MAKE) -C $(CUDD_DIR) -j$(NPROC)
	@test -f $@ || { echo "CUDD build finished but $@ is missing"; exit 1; }
	@echo "CUDD ready: $@"

# ---------------------------------------------------------------------------

config:
	@echo "CXX         $(CXX) ($(shell $(CXX) -dumpversion 2>/dev/null))"
	@echo "CXXFLAGS    $(CXXFLAGS)"
	@echo "BUILD_DIR   $(BUILD_DIR)"
	@echo "CUDD_INC    $(CUDD_INC)"
	@echo "CUDD_LIB    $(CUDD_LIB) [$(if $(wildcard $(CUDD_LIB)),found,MISSING)]"
	@echo "sources     $(words $(SRCS)) in src/"
	@echo "tests       $(words $(TEST_SRCS)) in test/"
	@echo "examples    $(words $(EX_SRCS)) in examples/"

install: $(LIB)
	@mkdir -p $(DESTDIR)$(PREFIX)/include/spbdd $(DESTDIR)$(PREFIX)/lib
	cp include/spbdd/*.hpp $(DESTDIR)$(PREFIX)/include/spbdd/
	cp $(LIB) $(DESTDIR)$(PREFIX)/lib/
	@echo "installed into $(DESTDIR)$(PREFIX)"

uninstall:
	rm -rf $(DESTDIR)$(PREFIX)/include/spbdd
	rm -f  $(DESTDIR)$(PREFIX)/lib/libspbdd.a

clean:
	rm -rf $(BUILD_DIR)

# Also throws away the CUDD build (but not the checked-out source tree).
distclean: clean
	@test -f $(CUDD_DIR)/Makefile && $(MAKE) -C $(CUDD_DIR) distclean || true

help:
	@echo "make            build $(LIB)"
	@echo "make cudd       fetch and build CUDD into $(CUDD_DIR)/"
	@echo "make check      build and run every test in test/"
	@echo "make examples   build everything in examples/"
	@echo "make config     show the paths and flags in use"
	@echo "make install    copy headers and library into PREFIX (=$(PREFIX))"
	@echo "make clean      remove $(BUILD_DIR)"
	@echo "make distclean  clean + undo the CUDD build"

-include $(OBJS:.o=.d)
