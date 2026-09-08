# ===========================================================================
#  spbdd -- a BDD-based library for sets of symplectic Pauli operators.
#
#  Layout
#    include/spbdd/   public headers   (installed)
#    src/             implementation   -> build/libspbdd.a
#    test/            one binary per .cpp, run by `make check`
#    examples/        one binary per .cpp, uses the public API only
#    buddy/           BuDDy source tree, built in place (not part of this repo)
#
#  Quick start
#    make buddy       fetch + build BuDDy into buddy/         (once)
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

# --- BuDDy -----------------------------------------------------------------
# Built in place under buddy/ -- no `make install`, nothing lands in
# /usr/local. BuDDy is a small set of .c files, so this compiles them directly
# and needs no autotools at all; the five-macro config.h its kernel.c expects
# is generated below.
BUDDY_REPO ?= https://github.com/utwente-fmt/buddy.git
BUDDY_DIR  ?= buddy
BUDDY_INC  ?= $(BUDDY_DIR)/src
BUDDY_GEN  ?= $(BUILD_DIR)/buddy
BUDDY_LIB  ?= $(BUILD_DIR)/libbdd.a

# --- sources ---------------------------------------------------------------
SRCS := $(wildcard src/*.cpp)
OBJS := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
LIB  := $(BUILD_DIR)/libspbdd.a

TEST_SRCS := $(wildcard test/*.cpp)
TEST_BINS := $(patsubst test/%.cpp,$(BUILD_DIR)/%,$(TEST_SRCS))
EX_SRCS   := $(wildcard examples/*.cpp)
EX_BINS   := $(patsubst examples/%.cpp,$(BUILD_DIR)/%,$(EX_SRCS))

INCLUDES := -Iinclude -I$(BUDDY_INC) -I$(BUDDY_GEN)
LDLIBS   := $(LIB) -lm
DEPFLAGS  = -MMD -MP

ifeq ($(wildcard $(BUDDY_LIB)),)
$(warning BuDDy not found at $(BUDDY_LIB) -- run `make buddy` first)
endif

# ===========================================================================

.PHONY: all lib buddy deps check examples config install uninstall clean distclean help
.DEFAULT_GOAL := all

all: lib

lib: $(LIB)

# libspbdd.a bundles BuDDy's objects as well, so a program that uses this
# library links one archive and needs no BuDDy on its include or library path.
$(LIB): $(OBJS) $(BUDDY_LIB)
	@mkdir -p $(dir $@) $(BUILD_DIR)/buddy_obj
	@rm -f $@ $(BUILD_DIR)/buddy_obj/*.o
	cd $(BUILD_DIR)/buddy_obj && $(AR) x $(abspath $(BUDDY_LIB))
	$(AR) rcs $@ $(OBJS) $(BUILD_DIR)/buddy_obj/*.o
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
buddy deps: $(BUDDY_LIB)

# BuDDy has an autotools build, but its kernel is eleven .c files that compile
# straight through -- so skip configure entirely and generate the only thing
# they need from it, which is five macros.
$(BUDDY_LIB):
	@command -v git >/dev/null 2>&1 || { \
	    echo "missing build tool: git"; \
	    echo "  sudo apt install -y build-essential git"; \
	    exit 1; }
	@test -d $(BUDDY_DIR) || git clone --depth 1 $(BUDDY_REPO) $(BUDDY_DIR)
	@mkdir -p $(BUDDY_GEN) $(BUILD_DIR)/buddy_build
	@printf '#define PACKAGE "buddy"\n#define PACKAGE_VERSION "2.4"\n#define VERSION "2.4"\n#define MAJOR_VERSION 2\n#define MINOR_VERSION 4\n' > $(BUDDY_GEN)/config.h
	@for f in $(BUDDY_DIR)/src/*.c; do \
	    $(CC) $(OPT) -I$(BUDDY_INC) -I$(BUDDY_GEN) -c $$f \
	        -o $(BUILD_DIR)/buddy_build/$$(basename $$f .c).o || exit 1; \
	done
	@$(AR) rcs $@ $(BUILD_DIR)/buddy_build/*.o
	@echo "BuDDy ready: $@ ($$($(AR) t $@ | wc -l) objects)"

# ---------------------------------------------------------------------------

config:
	@echo "CXX         $(CXX) ($(shell $(CXX) -dumpversion 2>/dev/null))"
	@echo "CXXFLAGS    $(CXXFLAGS)"
	@echo "BUILD_DIR   $(BUILD_DIR)"
	@echo "BUDDY_INC   $(BUDDY_INC)"
	@echo "BUDDY_LIB   $(BUDDY_LIB) [$(if $(wildcard $(BUDDY_LIB)),found,MISSING)]"
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

# BuDDy's objects live under build/, so `clean` already removed them.
distclean: clean
	@echo "nothing further to clean; remove $(BUDDY_DIR)/ by hand to re-fetch"

help:
	@echo "make            build $(LIB)"
	@echo "make buddy      fetch and build BuDDy into $(BUDDY_DIR)/"
	@echo "make check      build and run every test in test/"
	@echo "make examples   build everything in examples/"
	@echo "make config     show the paths and flags in use"
	@echo "make install    copy headers and library into PREFIX (=$(PREFIX))"
	@echo "make clean      remove $(BUILD_DIR)"
	@echo "make distclean  clean"

-include $(OBJS:.o=.d)
