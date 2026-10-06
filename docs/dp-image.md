# Deciding inequivalent pairs by a dynamic program over the error-set diagram

*Design note for `StabilizerCode::find_inequivalent_pair` (`PairMethod::Dp`).*

This note defines the problem, states the method, proves it correct, analyses its cost, and compares it with the
change-of-basis ("compose") implementation it replaces as the default. Section 8 gives measurements; all numbers in
it can be regenerated with `examples/pair_benchmark.cpp`.

---

## 1. Problem

Let $S \le \mathcal P_n$ be a stabilizer group with $r$ independent generators $g_1,\dots,g_r$, so the code has
$k = n - r$ logical qubits. Work with phase-free Paulis in symplectic form: $e = (x \mid z) \in \mathbb F_2^{2n}$, with
SPBDD's variable convention $x_q = $ variable $2q$, $z_q = $ variable $2q+1$. Write
$\langle e, f\rangle = x_e\!\cdot\! z_f + z_e\!\cdot\! x_f \pmod 2$ for the symplectic form.

Fix logical representatives $\bar X_1,\dots,\bar X_k,\bar Z_1,\dots,\bar Z_k \in N(S)$ with
$\langle \bar X_i,\bar Z_j\rangle = \delta_{ij}$ and $\langle \bar X_i,\bar X_j\rangle = \langle \bar Z_i,\bar Z_j\rangle = 0$.
Define two $\mathbb F_2$-linear maps

$$
\sigma(e) = \big(\langle e, g_i\rangle\big)_{i=1}^{r} \in \mathbb F_2^{r},
\qquad
\lambda(e) = \big(\langle e,\bar Z_j\rangle_{j=1}^{k},\ \langle e,\bar X_j\rangle_{j=1}^{k}\big) \in \mathbb F_2^{2k},
$$

the **syndrome** and the **logical signature**. In the code these are `syndrome()` and `logical_signature()`; the
concatenation $\pi = (\sigma,\lambda) : \mathbb F_2^{2n} \to \mathbb F_2^{m}$, $m = r + 2k = n + k$, is the projection of
the change of basis onto the coordinates that matter.

**Problem (inequivalent pair).** Given a set $\mathcal E \subseteq \mathbb F_2^{2n}$ of errors, represented as a
BDD over the $2n$ variables, decide whether there are $e_1, e_2 \in \mathcal E$ with $e_1 + e_2 \in N(S)\setminus S$
(and return such a pair).

**Lemma 1 (pairwise criterion).** For $e_1, e_2 \in \mathbb F_2^{2n}$:
$e_1 + e_2 \in N(S)\setminus S \iff \sigma(e_1)=\sigma(e_2)\ \wedge\ \lambda(e_1)\neq\lambda(e_2)$.

*Proof.* $N(S) = \{ f : \sigma(f) = 0\}$ and, for $f \in N(S)$, $f \in S \iff \lambda(f) = 0$ (an element of $N(S)$
is $\sum a_i g_i + \sum (c_j \bar X_j + f_j \bar Z_j)$, and its $\bar X,\bar Z$ coefficients are read off by pairing with the
partners $\bar Z_j,\bar X_j$). Linearity of $\sigma$ and $\lambda$ turns $\sigma(e_1+e_2)=0$ and $\lambda(e_1+e_2)\neq 0$
into the stated condition. $\square$

Hence, with the **relation**

$$
G(\mathcal E) \;=\; \pi(\mathcal E) \;=\; \{(\sigma(e),\lambda(e)) : e \in \mathcal E\} \subseteq \mathbb F_2^{r}\times\mathbb F_2^{2k},
$$

an inequivalent pair exists iff some syndrome $s$ occurs in $G$ with at least two different signatures.

**Lemma 2 (2k single-bit tests).** Such an $s$ exists iff for some signature bit $v \in\{1,\dots,2k\}$
$$
\exists \lambda.\,(G \wedge \lambda_v = 0)\ \ \wedge\ \ \exists \lambda.\,(G \wedge \lambda_v = 1)
$$
is satisfiable (as a function of $s$). *Proof.* Two different signatures differ in some bit $v$; conversely a syndrome
admitting both values of one bit has two signatures. $\square$

So the whole decision reduces to **building $G(\mathcal E)$** as a BDD over $r+2k$ variables; the $2k$ tests of
Lemma 2 are shared by every method (`find_inequivalent_pair`, "multi" loop). Everything below is about the cost of
that one step.

---

## 2. The compose method (reference)

Complete $\{g_i,\bar X_j,\bar Z_j\}$ to a symplectic basis with destabilizers $d_i$, and let $T$ be the invertible
$2n\times 2n$ matrix sending $e$ to its coordinates $(a\mid b\mid c\mid f)$ with $b=\sigma(e)$, $(c,f)=\lambda(e)$ and
$a$ the destabilizer-paired coefficients. The characteristic function of $T(\mathcal E)$ is
$\chi_{\mathcal E}\circ T^{-1}$, so one substitutes, **simultaneously for all $2n$ variables**, the variable $y_v$ of
$\chi_{\mathcal E}$ by the XOR of the new variables selected by row $v$ of $T^{-1}$ (`to_code_coordinates`), and then
existentially quantifies the $a$-block (`relation_by_compose`).

The cost of the substitution is governed by the size of its *result*, which is not controlled by the size of the
input. Example: $f(x,y)=\bigwedge_i \neg y_i$ has $n+2$ nodes; the substitution $y_i := x_i\oplus y_i$ turns it into
the equality function $\bigwedge_i (x_i\leftrightarrow y_i)$, which needs $\ge 2^n$ nodes in the order
$x_1..x_n,y_1..y_n$. A change of basis whose rows have many ones does this on a $2n$-variable scale, and the
stabilizer block is only discarded afterwards.

---

## 3. The DP method

### 3.1 Idea

Only the *image* $\pi(\mathcal E)$ is needed, and $\pi$ is linear. Write $e = \sum_{v=0}^{2n-1} e_v\,u_v$ with $u_v$ the
unit vectors (the single-qubit Paulis $X_q$ for $v=2q$, $Z_q$ for $v=2q+1$) and put

$$
w_v \;=\; \pi(u_v) \in \mathbb F_2^{m}\qquad(\text{the "flip vector" of variable } v),
\qquad\text{so}\qquad \pi(e) = \bigoplus_{v} e_v\, w_v .
$$

A BDD for $\mathcal E$ enumerates the assignments $e$ by Shannon decomposition on the variables in level order. The image
of a node is therefore built from the images of its two children by *adding one constant vector*, never by
re-expressing a variable as a XOR of other variables.

### 3.2 Definitions

Let the levels of the manager be $0,\dots,N-1$ ($N=2n$), $\mathrm{var}(\ell)$ the variable at level $\ell$, and for a
Boolean function $F$ over $\mathbb F_2^m$ and a vector $w$ define the **shift**

$$
(w \oplus F)(y) \;=\; F(y \oplus w),\qquad\text{i.e.}\qquad w\oplus F=\{\,w\oplus t : t\in F\,\}.
$$

For a node $u$ of the diagram of $\mathcal E$ and a level $\ell \le N$, let $\mathcal E_{u,\ell}$ be the set of
assignments of the variables at levels $\ell,\dots,N-1$ that $u$ accepts when the variables above level $\ell$ are
ignored (a variable skipped by the diagram is unconstrained). The **suffix image** is

$$
I(u,\ell) \;=\; \Big\{ \bigoplus_{\ell'\ge \ell} e_{\mathrm{var}(\ell')}\, w_{\mathrm{var}(\ell')} \;:\; e\in\mathcal E_{u,\ell} \Big\}\subseteq \mathbb F_2^m .
$$

Then $G(\mathcal E) = I(\text{root},0)$.

### 3.3 Recurrence

With $v=\mathrm{var}(\ell)$:

| case | $I(u,\ell)$ |
|---|---|
| $u=\bot$ | $\emptyset$ |
| $\ell=N$ ($u=\top$) | $\{0\}$ |
| $u$ is a node at level $\ell$ | $I(\mathrm{low}(u),\ell{+}1)\ \cup\ \big(w_v \oplus I(\mathrm{high}(u),\ell{+}1)\big)$ |
| otherwise (level $\ell$ is skipped, or $u=\top$ with $\ell<N$) | $I(u,\ell{+}1)\ \cup\ \big(w_v \oplus I(u,\ell{+}1)\big)$ |

**Theorem 3 (correctness).** The recurrence computes $I(u,\ell)$ as defined in 3.2. In particular
$I(\text{root},0)=\pi(\mathcal E)=G(\mathcal E)$.

*Proof.* Induction on $N-\ell$. For $\ell=N$ the only suffix assignment is empty and its image is $\{0\}$ if $u$
accepts it ($u=\top$) and $\emptyset$ otherwise. For $\ell<N$, split $\mathcal E_{u,\ell}$ by the value $b$ of the variable
at level $\ell$. If $u$ tests that variable, the assignments with $b=0$ are exactly $\mathcal E_{\mathrm{low}(u),\ell+1}$
and contribute images $I(\mathrm{low}(u),\ell+1)$; those with $b=1$ are $\mathcal E_{\mathrm{high}(u),\ell+1}$ and
contribute $w_v\oplus(\cdot)$ of their images, by linearity of $\bigoplus$. If $u$ does not test it, both values of $b$ leave
the remaining constraint $\mathcal E_{u,\ell+1}$ unchanged, giving the last row. $\square$

### 3.4 The shift is a relabelling

$w\oplus F$ is obtained from the diagram of $F$ by exchanging the two children of every node labelled by a variable in
$\mathrm{supp}(w)$. This is a bijection on nodes that preserves reducedness (two nodes are equal iff their images are,
and a node is redundant iff its image is), hence **$|w\oplus F| = |F|$ exactly**. The implementation expresses it as a
substitution of each flipped variable by its own negation (`Bdd::compose` with literals $\neg y_t$), never by an XOR of
several variables.

### 3.5 Algorithm

```
image(u, ℓ):                                  -- memo keyed on (u, ℓ)
    if u = ⊥:          return ∅
    if ℓ = N:          return {0}               -- the cube y = 0 on the m target variables
    if (u,ℓ) in memo:  return memo[(u,ℓ)]
    v ← var(ℓ)
    if u is a node and level(top(u)) = ℓ:
        R ← image(low(u), ℓ+1)  ∪  shift(image(high(u), ℓ+1), w_v)
    else:                                       -- skipped level
        H ← image(u, ℓ+1)
        R ← H ∪ shift(H, w_v)                   -- = H when w_v = 0
    memo[(u,ℓ)] ← R;  return R

find_inequivalent_pair(E):
    G ← image(root(E), 0)                       -- over the b, c, f variables
    run the 2k single-bit tests of Lemma 2 on G
    if none fires: return "no pair"
    otherwise pick an offending syndrome s and two signatures λ1 ≠ λ2 in the fibre G∧(σ=s),
    and return any members of  E ∩ with_syndrome(s) ∩ with_logical_signature(λ_i)
```

Notes on the implementation (`src/stabilizercode.cpp`).

* **Precomputation.** The constructor stores, for each of the $2n$ variables, the list of target variables
  (`b_var`, `c_var`, `f_var`) on which $w_v$ is 1 (`flips_`): $2n\cdot m$ symplectic products, $O(n^3)$ bit operations —
  no more than the linear algebra the constructor already does. The destabilizer block is not stored: it affects only
  the $a$ coordinates, which the quotient by $S$ discards.
* **Target variables** are the same `b_var/c_var/f_var` numbers the compose route uses, so both routes return the
  *identical canonical diagram* $G$ (tested), and the multiplicity tests and witness extraction are shared.
* **Witnesses** are now read in the original coordinates, `errors & with_syndrome(s) & with_logical_signature(λ)`,
  because $G$ forgets which error produced each point. This holds for both methods.
* **Reordering.** The DP interleaves reading the levels of $\mathcal E$ with building other diagrams. A dynamic
  reorder in between would move the levels, so `ReorderPause` switches dynamic reordering off for the duration of the call
  and restores it (`Manager::dynamic_reordering()` was added to make that possible).
* The recursion depth is at most $2n+1$; memo entries hold diagrams and are released when the call returns.

---

## 4. Complexity

Notation: $|f|$ is the node count of the input diagram, $N=2n$, $m=r+2k$, $\mathrm{skip}(f)=\sum_{(u\to c)}
(\mathrm{level}(c)-\mathrm{level}(u)-1)$ over the edges of $f$ (levels skipped).

* **Number of memo entries.** Every $(u,\ell)$ reached is either $(u,\mathrm{level}(u))$ or lies on a chain of skipped
  levels entered through an edge into $u$. Hence $M \le |f| + \mathrm{skip}(f) \le |f|\,(1+N)$, and $M\approx|f|$ for
  diagrams that skip few levels. The terminal $\top$ contributes at most $N$ entries.
* **Work per entry.** At most one shift and one union. A shift returns a diagram of exactly the same size (3.4) and, as a relabelling, needs one pass per flipped
  variable in the ideal case; the BuDDy `bdd_veccompose` used for it is not guaranteed to be linear, so only "polynomial in
  the diagram size" is claimed. A union of $A$ and $B$ costs $O(|A|\,|B|)$ and returns at most that many nodes.
* **Total.**
  $$T_{\mathrm{DP}} \;=\; \sum_{(u,\ell)\in\text{memo}} O\big(|I(\mathrm{lo})|\cdot|I(\mathrm{hi})|\ +\ |w_v|\,|I(\mathrm{hi})|\big),
  \qquad
  S_{\mathrm{DP}} \;=\; O\Big(\sum_{(u,\ell)} |I(u,\ell)|\Big).$$
  Every $I(u,\ell)$ is a subset of $\mathbb F_2^m$, so each diagram has at most $2^{m}$ points, and $I(\mathrm{root},0)=G$.
  The bound is *output-sensitive in the suffix images*: it is small whenever the images of suffixes of $\mathcal E$ stay
  structured in $(\sigma,\lambda)$ space, and there is no a-priori polynomial bound — the decision problem is not claimed
  to be tractable in general.
* **Multiplicity tests** (shared): $2k$ rounds of two conjunctions and two quantifications over $2k$ variables, each
  $O(|G|)$ up to the apply cost; negligible next to building $G$ in every measurement below.
* **Compose route** for comparison: one simultaneous substitution of $N$ variables by XORs whose total length is the
  number of ones of $T^{-1}$ (up to $N^2$), followed by $\exists a$. Its running time is a function of the size of the
  substituted diagram on $2n$ variables, which can be exponentially larger than $|f|$ and than $|G|$ (Section 2).

Why the DP wins in practice, in one line: both methods end at the same $G$ on $r+2k=n+k$ variables, but the compose route
must pass through a $2n$-variable diagram coupled by XORs, while the DP only ever holds diagrams over the $n+k$ target
variables and moves between them by relabelling.

---

## 5. Relation to the per-location construction

For an error set built as a sum of small sets — "at most $t$ of these locations, each contributing one of a few Paulis" —
the image can also be built by layered sets, $G_j \leftarrow G_j \cup \bigcup_p (\pi(p)\oplus G_{j-1})$, with no input
diagram at all. On weight balls the node counts of that construction and of the DP coincide ($G$ is canonical), and the
DP needs no knowledge of how $\mathcal E$ was built: it applies to **any** `PauliSet`.

---

## 6. API changes

```cpp
enum class StabilizerCode::PairMethod { Dp, Compose };
bool has_inequivalent_pair (const PauliSet&, PairMethod = PairMethod::Dp) const;
std::optional<Pair> find_inequivalent_pair(const PauliSet&, PairMethod = PairMethod::Dp) const;
bool Manager::dynamic_reordering() const;     // new, read-only
```

The compose implementation is kept (`PairMethod::Compose`) as the reference and for benchmarking.

## 7. Testing

`test/stabilizercode_test.cpp`, section *"the DP and the compose route agree"*: 60 random unions of cosets for each of the
bit-flip, five-qubit and Steane codes (both safe and unsafe sets occur), verdicts compared between methods, and every
witness from either method verified independently (members of the set, equal syndromes, different signatures, product
in $N(S)\setminus S$); weight balls with known answers; behaviour with dynamic reordering switched on. All
pre-existing tests (brute force on all subsets, the product formulation $E{*}E\cap(N(S)\setminus S)$) pass unchanged.

## 8. Measurements

All runs: one cloud container (2 vCPU, 7 GB), `g++ -O2`, BuDDy node table $2^{22}$ nodes / $2^{20}$ cache, default
variable order, no reordering, each row one run (no repeats). "compose" is `PairMethod::Compose`, "DP" is
`PairMethod::Dp`; both return the same verdict in every row where both finished. Two jobs were sometimes running on the two
cores at once, which slows both (the $d{=}11,t{=}4$ row took 112.8 s alone and 257 s with another job running), so read the
timings to a factor of about two. Table size matters as much: with BuDDy's smaller default table ($2^{20}$) compose on
$d{=}7,t{=}3$ did not finish in 600 s although it needs 111 s at $2^{22}$.

**(a) Weight balls** — every Pauli of weight $\le t$ on the rotated surface code (answer: no inequivalent pair, $d>2t$).

| $d$ | $t$ | $|\mathcal E|$ | diagram nodes | DP | compose |
|---|---|---|---|---|---|
| 5 | 2 | 2,776 | 139 | 0.011 s | 0.026 s |
| 7 | 2 | 1.07e4 | 283 | 0.037 s | 1.55 s |
| 7 | 3 | 5.08e5 | 369 | 0.176 s | 111.1 s |
| 9 | 3 | 2.33e6 | 625 | 0.916 s | not finished in 600 s |
| 9 | 4 | 1.37e8 | 771 | 6.16 s | not run |
| 11 | 4 | 6.96e8 | 1,171 | 112.8 s | not run |
| 13 | 4 | 2.7e9 | not measured | not finished in 900 s | not run |

**(b) Unstructured sets** — union of $K$ random cosets, each of $G$ random generators of weight $\le W$, on the same code
(`pair_benchmark cosets`; fixed seed). Answer: no pair in all four.

| $d$ | $K,G,W$ | $|\mathcal E|$ | diagram nodes | DP | compose |
|---|---|---|---|---|---|
| 5 | 20, 6, 4 | 1,280 | 7,265 | 0.021 s | 0.289 s |
| 7 | 50, 8, 5 | 1.28e4 | 113,772 | 0.522 s | 65.2 s |
| 9 | 100, 8, 5 | 2.56e4 | 387,759 | 5.07 s | not finished in 600 s |
| 11 | 300, 8, 5 | 7.67e4 | 1,569,997 | 257 s (two jobs running) | not run |

**(c) Error sets of a circuit** — all faults (at most $t$ faulty locations; input data errors counted as faults) of one
round of flag-qubit syndrome extraction, restricted to the all-flags-zero outcomes, with every measurement record kept as a
Pauli component on the measured qubit (so the code handed to `StabilizerCode` is $S\otimes I$ plus one single-qubit
generator per measured qubit, $N=n+\#\text{ancilla}+\#\text{flag}$ qubits). These were produced by a separate
experimental harness around the library, not by code in this repository; the fault set was built with `cx` / `fault_inject`
and `|`. Time is that of `find_inequivalent_pair` only.

| circuit | $N$ | $t$ | set nodes | DP | compose | explicit hash of every fault combination |
|---|---|---|---|---|---|---|
| Steane, flagged | 19 | 2 | 7,405 | 0.012 s | 0.125 s | 0.002 s |
| Golay [[23,1,7]], flagged | 67 | 1 | 11,591 | 0.022 s | 550 s | 0.004 s |
| Surface-5, flagged | 65 | 2 | 322,470 | 0.708 s | 56.3 s | 0.023 s |
| Golay, flagged | 67 | 2 | ~4e6 (no reordering) | 16.4 s | not completed | 0.087 s |

**What the numbers say.**

* On every set where compose finishes, DP is faster, from about 2x (smallest sets) to about 25,000x (Golay, $t{=}1$:
  550 s vs 0.022 s). On the $d{=}7,t{=}3$ ball it is 630x, on the flagged surface-5 set 80x. Sets that compose could not finish in 600 s are solved by DP in
  seconds ($d{=}9$ ball 0.9 s, $d{=}9$ cosets 5 s).
* DP does not make the diagram method competitive with a plain hash of the fault combinations when the set is small
  enough to list (rows of (c), last column). Its advantage over explicit enumeration only appears for sets with $10^8$ or more
  elements (the large weight balls), and there only moderately: an explicit multi-pass hash needed 252 s on the
  $d{=}11,t{=}4$ ball, against 113 s here.
* The in-library DP over the *diagram of the ball* is slower on the largest ball than a layered construction that builds
  $G$ directly from the per-qubit images (34.6 s for $d{=}11,t{=}4$ in an experimental harness, measured with another job
  running). That construction needs to know that the set is "at most $t$ of these locations"; the DP does not. The
  difference was not analysed.
* The multiplicity tests are identical in both methods; the timings above include them, and the large differences come from
  how $G$ is obtained (for compose, the separately timed substitution accounted for 99% or more of the time on the $d{=}7,t{=}3$
  ball and on Golay $t{=}1$).


## 9. Limitations and open points

* No polynomial bound on the suffix images; a set whose suffix images are unstructured in $(\sigma,\lambda)$ space can still
  be expensive. The measurements show where it is not.
* The DP follows the level order of the *input* diagram. How much a better input order helps was not studied.
* Reordering is paused, not exploited, inside the DP.
* For error sets small enough to list, hashing $(\sigma,\lambda)$ per element is far cheaper than any diagram method (table (c)).
  The diagram methods are for sets too large to list.
