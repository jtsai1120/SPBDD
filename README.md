# SPBDD: A BDD-based C++ library for operations on sets of Pauli operators.

SPBDD stores sets of *phase-free* $n$-qubit Pauli operators in symplectic representation as binary decision diagrams and provides set operations, Clifford gate propagation, fault injection, nontrivial logical error pair detection, etc.

For example, finding the distance of a code:

```cpp
#include <spbdd/spbdd.hpp>
using namespace spbdd;

PauliSpace sp(7);                                  // 7 qubits

// The Steane [[7,1,3]] code.
std::vector<std::string> generators = {"IIIXXXX", "IXXIIXX", "XIXIXIX",
                                       "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"};

PauliSet S = sp.generated_by(generators);          // the stabilizer group
PauliSet N = sp.commuting_with_all(generators);    // its normalizer

printf("d = %d\n", (N - S).min_weight());          // the code distance: 3
```

Moreover, we can inject two-qubit gate faults, propagating the error set through a circuit, and asking whether a
decoder can still tell every pair of the resulting errors apart (i.e. nontrivial logical error pair detection):

```cpp
StabilizerCode code(sp, generators);

// Create a pauli set {IIIIIII} and add a two-qubit fault on qubit 0, 1
PauliSet errors = sp.identity().fault_inject({0, 1});

// Propagate the pauli set through the gates
errors = errors.cx(0, 2).h(3).cz(3, 5);

if (auto pair = code.find_inequivalent_pair(errors))
    printf("unsafe: %s and %s forms a nontrivial logical error pair\n",   
           pair->first.c_str(), pair->second.c_str()); // IIIIIII and XXXIIII
```

## Installation and Compilation

SPBDD is built on [CUDD](https://github.com/ivmai/cudd) and requires a C++17
compiler. On Ubuntu/Debian:

```bash
sudo apt install -y build-essential git autoconf automake libtool
```

Then, from the top of the repository:

```bash
make cudd      # clone and build CUDD into cudd/ (once)
make           # build build/libspbdd.a
make check     # build and run the test suite
```

To install the headers and the archive system-wide:

```bash
sudo make install                   # into /usr/local by default
make install PREFIX=$HOME/.local    # or somewhere else
```

## Usage

Include the single header and link the single archive:

```cpp
#include <spbdd/spbdd.hpp>
```

```bash
g++ -std=c++17 -Iinclude my_program.cpp build/libspbdd.a -lm -o my_program
```

After `make install` the include path is no longer needed and the library is
found by name:

```bash
g++ -std=c++17 my_program.cpp -lspbdd -lm -o my_program
```

A worked example covering the whole library — code data, set algebra, gate
propagation, fault injection, measurement, and the decoder-safety check
(logical error pair detection) — is in
[`examples/demo.cpp`](examples/demo.cpp):

```bash
make examples && ./build/demo
```



## API Reference

Everything lives in the single header `<spbdd/spbdd.hpp>` and in the namespace `spbdd::` .

To use the tool, we always start from creating a `PauliSpace` instance, which is the workspace and factory of `PauliSet`. 

```cpp
spbdd::PauliSpace(int n_qubits);
```

Then, we create `PauliSet` instance by `PauliSpace`'s public functions and do set operations by `PauliSet`'s public functions.

Note that sets built from the same `PauliSpace` may be operated with each other, while sets from different `PauliSpace` may not.

### `PauliSpace` — building set(s)



| Returns | Member | Arguments | Description |
|---|---|---|---|
| `int` | `n_qubits` | | Number of qubits |
| `void` | `grow_to` | `int n_qubits` | Add qubits; the space can only grow |

**Basic sets**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `empty` | | `{}` |
| `PauliSet` | `all` | | All 4^n operators |
| `PauliSet` | `identity` | | `{ I...I }` |
| `PauliSet` | `from` | `const std::string &s` | The single operator `s`, e.g. `"IXYZ"` |
| `PauliSet` | `from_list` | `const std::vector<std::string> &list` | Exactly the listed operators |

**Subgroups**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `generated_by` | `const std::vector<std::string> &generators` | The group the generators span, 2^rank elements |
| `PauliSet` | `coset_of` | `const std::string &base, const std::vector<std::string> &generators` | `base` multiplied by that group |

**Support**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `supported_on` | `const std::vector<int> &qubits` | The identity outside these qubits, free on them |
| `PauliSet` | `identity_on` | `const std::vector<int> &qubits` | The identity on these qubits, free elsewhere |
| `PauliSet` | `pauli_at` | `int qubit, char p` | Exactly `p` on that qubit, free elsewhere |
| `PauliSet` | `non_identity_at` | `int qubit` | Anything but the identity on that qubit |
| `PauliSet` | `matching` | `const std::string &pattern` | A glob: `I X Y Z` pin a qubit, `*` `?` `.` leave it free |

**Weight**

The weight of an operator is the number of qubits it does not act on with the
identity. The overloads taking `qubits` count only those, leaving the rest
unconstrained.

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `weight_exactly` | `int w` | Operators of weight exactly `w` |
| `PauliSet` | `weight_at_most` | `int w` | Operators of weight at most `w` |
| `PauliSet` | `weight_between` | `int lo, int hi` | Operators of weight in `[lo, hi]` |
| `PauliSet` | `weight_exactly` | `int w, const std::vector<int> &qubits` | |
| `PauliSet` | `weight_at_most` | `int w, const std::vector<int> &qubits` | |
| `PauliSet` | `weight_between` | `int lo, int hi, const std::vector<int> &qubits` | |

**Commutation**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `commuting_with` | `const std::string &p` | Everything commuting with `p` |
| `PauliSet` | `anticommuting_with` | `const std::string &p` | Everything anticommuting with `p` |
| `PauliSet` | `commuting_with_all` | `const std::vector<std::string> &ps` | Everything commuting with all of them; for stabilizer generators, the normaliser |

### `PauliSet` — operation on set(s)

Operands of a binary operation must come from the same `PauliSpace`. Every
member returns a new set rather than modifying this one, apart from the compound
assignments.

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `PauliSpace` | `space` | | The universe this set lives in |
| `int` | `n_qubits` | | Number of qubits |

**Set algebra**

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `PauliSet` | `operator\|` | `const PauliSet &o` | Union |
| `PauliSet` | `operator&` | `const PauliSet &o` | Intersection |
| `PauliSet` | `operator-` | `const PauliSet &o` | Difference |
| `PauliSet` | `operator^` | `const PauliSet &o` | Symmetric difference |
| `PauliSet` | `operator~` | | Complement |
| `PauliSet &` | `operator\|=` `operator&=` `operator-=` `operator^=` | `const PauliSet &o` | The same, in place |
| `bool` | `operator==` `operator!=` | `const PauliSet &o` | Set equality |

Note that C++ gives `|`, `&` and `^` lower precedence than `==`, so a comparison
against a combination needs parentheses: `x == (a ^ b)`.

**Predicates**

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `bool` | `is_empty` | | |
| `bool` | `is_universe` | | Whether the set is all 4^n operators |
| `bool` | `contains` | `const std::string &p` | |
| `bool` | `subset_of` | `const PauliSet &o` | |
| `bool` | `superset_of` | `const PauliSet &o` | |
| `bool` | `disjoint_from` | `const PauliSet &o` | |

**Metrics**

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `double` | `size` | | Number of elements; a `double` because it reaches 4^n |
| `std::size_t` | `node_count` | | Number of nodes in the diagram, which is what the set costs |
| `int` | `min_weight` | | The smallest weight present, or `-1` for the empty set |

**Derived sets**

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `PauliSet` | `fault_inject` | `const std::vector<int> &qubits` | Inject a fault at these locations: multiply by every operator supported on them |
| `PauliSet` | `reset` | `const std::vector<int> &qubits` | Rewrite every element to be the identity on these qubits |
| `PauliSet` | `with_weight_at_most` | `int w` | The elements of weight at most `w` |
| `PauliSet` | `with_weight_exactly` | `int w` | The elements of weight exactly `w` |

**Clifford gates**

Every gate acts by conjugation, `E -> U E U*`. Because no phase is tracked, the
`X`, `Y` and `Z` gates are the identity map; they are provided so a circuit can
be replayed without special-casing them.

| Returns | Member | Arguments | Images of the generators |
|---|---|---|---|
| `PauliSet` | `x` `y` `z` | `int qubit` | unchanged |
| `PauliSet` | `h` | `int qubit` | `X -> Z`, `Z -> X` |
| `PauliSet` | `s` | `int qubit` | `X -> Y`, `Z -> Z`; `S` and `S*` agree here |
| `PauliSet` | `sx` | `int qubit` | `X -> X`, `Z -> Y` |
| `PauliSet` | `cx` | `int control, int target` | `X_c -> X_c X_t`, `Z_t -> Z_c Z_t` |
| `PauliSet` | `cy` | `int control, int target` | `X_c -> X_c Y_t`, `X_t -> Z_c X_t`, `Z_t -> Z_c Z_t` |
| `PauliSet` | `cz` | `int a, int b` | `X_a -> X_a Z_b`, `X_b -> Z_a X_b` |
| `PauliSet` | `swap` | `int a, int b` | The two qubits exchange |

**Measurement**

Measuring splits a set rather than changing it: an element either commutes with
the observable, and the outcome is the one the ideal circuit would give, or it
anticommutes and the outcome is flipped. An empty side is an outcome that cannot
occur. Measuring a stabilizer generator is one bit of syndrome extraction. A
destructive measurement that reuses the qubit is a split followed by `reset()`.

```cpp
struct MeasurementSplit {
    PauliSet unflipped;   // the outcome the ideal circuit would give
    PauliSet flipped;     // the other one
};
```

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `MeasurementSplit` | `measure` | `const std::string &observable` | Measure any Pauli observable |
| `MeasurementSplit` | `measure_z` | `int qubit` | Z-basis measurement |
| `MeasurementSplit` | `measure_x` | `int qubit` | X-basis measurement |
| `MeasurementSplit` | `measure_y` | `int qubit` | Y-basis measurement |

**Enumeration**

These are the only operations whose cost is the number of elements.

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `void` | `for_each` | `const std::function<bool(const std::string &)> &fn` | Call `fn` once per element; return `false` from it to stop. The order is unspecified |
| `std::vector<std::string>` | `to_strings` | `std::size_t max_elements = 4096` | Every element, sorted; throws `std::length_error` beyond the limit |
| `std::optional<std::string>` | `any_element` | | Any one element, or nothing if the set is empty |

### `StabilizerCode` — code data and decoder safety

Built once from a set of generators and reused across queries. The generators
need not be independent, but they must commute with each other.

```cpp
StabilizerCode(PauliSpace space, const std::vector<std::string> &generators);
```

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `PauliSpace` | `space` | | The universe |
| `int` | `n_qubits` | | *n* |
| `int` | `n_stabilizers` | | *r*, the rank of the generator set |
| `int` | `n_logical` | | *k = n - r* |
| `const std::vector<std::string> &` | `stabilizers` | | An independent set of generators, *r* of them |
| `const std::vector<std::string> &` | `destabilizers` | | Their partners, *r* of them |
| `const std::vector<std::string> &` | `logical_x` | | Logical `X` operators, *k* of them |
| `const std::vector<std::string> &` | `logical_z` | | Logical `Z` operators, *k* of them |
| `PauliSet` | `group` | | The stabilizer group *S* |
| `PauliSet` | `normalizer` | | *N(S)* |
| `int` | `distance` | | The weight of the lightest element of *N(S) \ S* |
| `std::vector<bool>` | `syndrome` | `const std::string &e` | *r* bits: which stabilizers the error anticommutes with |
| `std::vector<bool>` | `logical_signature` | `const std::string &e` | *2k* bits identifying the logical class |
| `PauliSet` | `with_syndrome` | `const std::vector<bool> &syndrome` | Every operator with that syndrome |
| `PauliSet` | `with_logical_signature` | `const std::vector<bool> &signature` | Every operator with that signature |
| `bool` | `has_inequivalent_pair` | `const PauliSet &errors` | Whether two elements multiply into a logical operator |
| `std::optional<Pair>` | `find_inequivalent_pair` | `const PauliSet &errors` | The same question, with a counterexample |

A pair reported by `find_inequivalent_pair` is indistinguishable to any decoder:
the two errors produce identical observations, yet correcting one leaves a
logical error in the other case.

```cpp
struct StabilizerCode::Pair {
    std::string       first, second;                    // two members of `errors`
    std::vector<bool> syndrome;                         // the one they share
    std::vector<bool> signature_first, signature_second;
};
```

## Citation


## Contact

Please [create an issue](https://github.com/jtsai1120/SPBDD/issues) for bug
reports and feature requests, or contact us directly at *jtsai1120@gmail.com*.
