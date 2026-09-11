# SPBDD: A BDD-based C++ library for operations on sets of Pauli operators.

SPBDD stores sets of *phase-free* $n$-qubit Pauli operators in symplectic representation as binary decision diagrams (currently built on [BuDDy](https://github.com/utwente-fmt/buddy)) and provides set operations, Clifford gate propagation, fault injection, nontrivial logical error pair detection, etc.

For example, given an initial Pauli error set, we can inject two-qubit gate faults, propagating the set through a circuit, and asking whether there exist two errors in the set such that their product is a non-trivial logical error (which is important when constructing a decoder):

```cpp
#include <iostream>
#include <spbdd/spbdd.hpp>

// Declare a 7-qubit Pauli Space
spbdd::PauliSpace sp(7);

// Create the identity pauli set {IIIIIII} then inject a two-qubit fault on qubit 0, 1
spbdd::PauliSet errors = sp.identity().fault_inject({0, 1});

// Propagate the pauli set through the gates
errors = errors.cx(0, 2).h(3).cz(3, 5);

// Declare the Steane [[7,1,3]] code
std::vector<std::string> generators = {
    "IIIXXXX", "IXXIIXX", "XIXIXIX",
    "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"
};
spbdd::StabilizerCode code(sp, generators);

// Find if there exist non-trivial logical error pairs (i.e. inequivalent pair under same syndrome)
if (auto pair = code.find_inequivalent_pair(errors))
    std::cout << pair->first << " and " << pair->second << " form a nontrivial logical error pair\n";
    // stdout: IIIIIII and XXXIIII form a nontrivial logical error pair
```

## Installation and Compilation

On Ubuntu/Debian, the following packages are required:

```bash
sudo apt install -y build-essential git
```

SPBDD is available for the *installed* version and the *local build* version.

For both versions, first build the library from the top of the repository:

```bash
make buddy     # clone and build BuDDy into buddy/
make           # build build/libspbdd.a
make check     # build and run the test suite
```

Then, for the *installed* version:

```bash
sudo make install                   # into /usr/local by default
make install PREFIX=$HOME/.local    # or somewhere else
```

## Usage

Include the header in your program.

```cpp
#include <spbdd/spbdd.hpp>
```

For the *local build* version, copy the whole `SPBDD/` folder into your workspace folder, then compile:

```bash
g++ -std=c++17 -ISPBDD/include my_program.cpp SPBDD/build/libspbdd.a -lm -o my_program
```

For the *installed* version, compile directly:

```bash
g++ -std=c++17 my_program.cpp -lspbdd -lm -o my_program
```

(Note: `-lm` links `libm`, which BuDDy uses for pow/log)

A worked example covering most of the library functions is in [`examples/demo.cpp`](examples/demo.cpp):

```bash
make examples && ./build/demo
```



## API Reference

Everything lives in the namespace `spbdd::`. 

Start from creating a `PauliSpace` instance, which owns the $n$-qubit workspace and is utilized to construct basic `PauliSet` instances.

```cpp
spbdd::PauliSpace(int n_qubits);
```

> **Note: Only one `PauliSpace` may be alive at a time** — BuDDy keeps its variable
table in process-wide globals, so a second one throws `std::runtime_error`.

---

### `PauliSpace`

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `int` | `n_qubits` | | # of qubits |
| `void` | `grow_to` | `int n_qubits` | Add qubits (# of qubits can only grow) |



#### **Basic sets**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `empty` | | {} |
| `PauliSet` | `all` | | {$I,X,Y,Z$}$^{\otimes n}$ |
| `PauliSet` | `identity` | | {$I$}$^{\otimes n}$ |
| `PauliSet` | `from` | `const std::string &s` | A single operator (e.g. `s="IXXYZ"` $\rightarrow$ {$IXXYZ$}) |
| `PauliSet` | `from` | `const std::vector<std::string> &list` | Every listed operator (e.g. `list={"IXXYZ","ZZZZZ"}` $\rightarrow$ {$IXXYZ, ZZZZZ$}) |

#### **Groups**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `generated_by` | `const std::vector<std::string> &generators` | The group the generators span |
| `PauliSet` | `coset_of` | `const std::string &base, const std::vector<std::string> &generators` | `base` multiplied by that group |

#### **Support**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `supported_on` | `const std::vector<int> &qubits` | {$I,X,Y,Z$} on these qubits, while $I$ on the others |
| `PauliSet` | `identity_on` | `const std::vector<int> &qubits` | $I$ on these qubits, while {$I,X,Y,Z$} on the others |
| `PauliSet` | `pauli_at` | `int qubit, char p` | `p` on that qubit, while {$I,X,Y,Z$} on the others |
| `PauliSet` | `non_identity_at` | `int qubit` | {$X,Y,Z$}  on that qubit, while {$I,X,Y,Z$} on the others |
| `PauliSet` | `matching` | `const std::string &pattern` | A specific set pattern, using `*` to represent {$I,X,Y,Z$} |

#### **Weight**

The weight of an operator is the number of qubits with {$X,Y,Z$}.

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `weight_exactly` | `int w [, const std::vector<int> &qubits]` | Weight exactly `w` |
| `PauliSet` | `weight_at_most` | `int w [, const std::vector<int> &qubits]` | Weight at most `w` |
| `PauliSet` | `weight_between` | `int lo, int hi [, const std::vector<int> &qubits]` | Weight in `[lo, hi]` |

> **Note:** `qubits` is an optional argument for constraint on certain qubits while counting weight.  

#### **Commutation**

| Returns | Member | Arguments | Resulting set |
|---|---|---|---|
| `PauliSet` | `commuting_with` | `const std::string &p` | Everything commuting with `p` |
| `PauliSet` | `anticommuting_with` | `const std::string &p` | Everything anticommuting with `p` |
| `PauliSet` | `commuting_with_all` | `const std::vector<std::string> &ps` | Everything commuting with all of them; for stabilizer generators, the normaliser |

---

### `PauliSet`

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `PauliSpace` | `space` | | The PauliSpace it lives in |
| `int` | `n_qubits` | | # of qubits |

#### **Set algebra**

> **Hint:** one can instead think of the boolean algebra on the characteristic functions instead of the set algebra.

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `PauliSet` | `operator\|` | `const PauliSet &o` | Union |
| `PauliSet` | `operator&` | `const PauliSet &o` | Intersection |
| `PauliSet` | `operator-` | `const PauliSet &o` | Difference |
| `PauliSet` | `operator^` | `const PauliSet &o` | Symmetric difference |
| `PauliSet` | `operator~` | | Complement |
| `PauliSet&` | `operator\|=`, `operator&=`, `operator-=`, `operator^=` | `const PauliSet &o` | In place ver. |
| `bool` | `operator==`, `operator!=` | `const PauliSet &o` | Set equality |

#### **Predicates**

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `bool` | `is_empty` | | Whether the set has no elements |
| `bool` | `is_universe` | | Whether the set is all $4^n$ operators |
| `bool` | `contains` | `const std::string &p` | Whether `p` is an element |
| `bool` | `subset_of` | `const PauliSet &o` | Whether every element is also in `o` |
| `bool` | `superset_of` | `const PauliSet &o` | Whether every element of `o` is also here |
| `bool` | `disjoint_from` | `const PauliSet &o` | Whether the two share no element |

#### **Metrics**

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `double` | `size` | | Number of elements; a `double` because it reaches $4^n$ |
| `std::size_t` | `node_count` | | Number of nodes in the diagram, which represents memery costs |
| `int` | `min_weight` | | The smallest weight present, or `-1` for the empty set |

#### **Derived sets**

| Returns | Member | Arguments | Description |
|---|---|---|---|
| `PauliSet` | `fault_inject` | `const std::vector<int> &qubits` | Inject a fault at these qubits |
| `PauliSet` | `reset` | `const std::vector<int> &qubits` | Rewrite every element to be the identity on these qubits |
| `PauliSet` | `with_weight_at_most` | `int w` | Filter the elements of weight at most `w` |
| `PauliSet` | `with_weight_exactly` | `int w` | Filter the elements of weight exactly `w` |

#### **Clifford gates**

Every gate $U$ acts on the Pauli operator $E$ by conjugation: $E \mapsto U E U^\dagger $.

| Returns | Member | Arguments | Images of the Pauli generators |
|---|---|---|---|
| `PauliSet` | `x` `y` `z` | `int qubit` | unchanged |
| `PauliSet` | `h` | `int qubit` | $X \mapsto Z,\ Z \mapsto X$ |
| `PauliSet` | `s`, `sdg` | `int qubit` | $X \mapsto Y,\ Z \mapsto Z$ |
| `PauliSet` | `sx`, `sxdg` | `int qubit` | $X \mapsto X,\ Z \mapsto Y$ |
| `PauliSet` | `cx` | `int control, int target` | $X_c \mapsto X_c X_t,\ Z_t \mapsto Z_c Z_t$ |
| `PauliSet` | `cy` | `int control, int target` | $X_c \mapsto X_c Y_t, X_t \mapsto Z_c X_t,\ Z_t \mapsto Z_c Z_t$ |
| `PauliSet` | `cz` | `int a, int b` | $X_a \mapsto X_a Z_b,\ X_b \mapsto Z_a X_b$ |
| `PauliSet` | `swap` | `int a, int b` | The two qubits exchange |

#### **Measurement**

Measuring splits a set rather than changing it: an element either commutes with
the observable and keeps the outcome the ideal circuit would give, or
anticommutes and gets the flipped one. An empty side is an outcome that cannot
occur. Measuring a stabilizer generator is one bit of syndrome extraction, and a
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

#### **Enumeration**

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
| `PauliSpace` | `space` | | The universe this code lives in |
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
| `std::optional<StabilizerCode::Pair>` | `find_inequivalent_pair` | `const PauliSet &errors` | The same question, with a counterexample |

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
