# Reproducing the finite computations

Archived with release v0.3.0 (see `../README.md`): the paths below are those of
the v0.2.0 layout, in which this directory was `computations/` at the
repository root; the certificate verifier is now `go/cmd/d6-certificate`.

Run from the repository root of release v0.2.0 with Python 3.10 or later:

```sh
python3 computations/run_search.py
python3 computations/verify_certificate.py
```

The code uses only the standard library. Equivalently, use
`uv run --no-project computations/run_search.py` and
`uv run --no-project computations/verify_certificate.py` if `uv` is available.
The first command regenerates `results/computation.json` and
`results/d6-recurrence-certificate.json`; the second verifies the latter
without calling the KL evaluator. The first command also runs this verifier.

## Conventions and algorithm

Elements are signed one-line permutations, and words multiply on the right.
Generator labels start at zero. In type D, generator 0 sends the first two
entries `(a,b)` to `(-b,-a)`, while generator `i>0` swaps the entries at
one-based positions `i,i+1`. Thus code labels `0,1,2,...` correspond to
Gern's `s_1,s_2,s_3,...`. In B, generator 0 negates the first entry, and
the other generators have the same swap convention. In A, generator i
swaps zero-based positions i and i+1. Polynomial arrays list coefficients
in ascending powers of q. All calculations use exact Python integers.

`coxeter.py` enumerates each group by breadth-first multiplication by its
generators. Reduced words, lengths, and the subword criterion give the
Bruhat lower ideals. Fully commutative status is computed recursively:
an element is non-FC when a right weak predecessor is non-FC or a reduced
expression ends in a full noncommuting braid. Every non-FC element receives
a concrete reduced-word witness containing such a braid.

The equal-parameter KL recurrence is Proposition 1.2.6 of Gern's thesis,
transported from left to right by inversion. If `ws<w` and `c=1` exactly
when `xs<x`, it is

```text
P(x,w) = q^(1-c) P(xs,ws) + q^c P(x,ws)
         - sum_{x<=z<ws, zs<z} mu(z,ws) q^((l(w)-l(z))/2) P(x,z).
```

`sparse_kl.py` supplies a second implementation. It avoids group
enumeration, computes lengths by signed inversion counts, and compares
Bruhat order by the lifting property. It is cross-checked against every
ordered pair in the five finite groups. Its KL polynomials are checked
against every comparable pair in the enumerated implementation.

The interval search uses deterministic directed color refinement to
reject impossible matches, then exact Hasse-graph backtracking. A separate
certificate check confirms the bijection and **every pairwise order
relation**, all subinterval KL polynomial equalities, all choices of right
descent in their recurrences, and their R-polynomial reciprocity identities.

## Scope and counts

The five exhaustive groups A3, A4, B3, D4, and B4 contain respectively
24, 120, 48, 192, and 384 elements. There are respectively 213, 3,781,
847, 9,817, and 40,249 comparable pairs, including equal endpoints:
54,907 in total. The independent Bruhat comparison checks cover all
201,600 ordered pairs across these groups.

All polynomial constant terms, degree bounds, positivity, invariance
under inversion, and `P(x,w0)=1` are checked. The fully commutative element
counts are 14, 42, 24, 48, and 83. The observed mu-values at strict FC-bottom
pairs are 0 or 1, as recorded individually by group in the JSON.

Six explicit interval certificates are saved: one B3 example, three B4
examples, and two A4 examples. Every target bottom is non-FC. Five
polynomials are `1+q`, at interval ranks 3 or 5. For this polynomial,
rank 3 has mu=1 and rank 5 has mu=0, despite its leading coefficient 1.
The sixth example is a rank-5 D4-to-B4 transfer with `P=1+2q+q^2`,
mu=1, and rank vector `(1,4,10,12,6,1)`. Its source top uses all of the
D4 fork. Its endpoints are source `(-2,-1,3,4)` to `(-4,-3,-2,-1)`,
target `(2,4,3,-1)` to `(2,3,-1,-4)`.
The B certificates include pairs of cover reflections whose product has
order 4. This excludes containment of the whole interval in a coset of
any type-A parabolic subgroup: two reflections in a type-A group have
product order 1, 2, or 3. Both left- and right-coset obstructions are saved.

These examples illustrate transfer and the failure of FC status to be
an abstract interval property. Their small ranks mean they do **not**
require the supplied unrestricted invariance theorem, and the finite
search proves no general transfer theorem or whole-type coverage.

## D6 base calculation and certificate

The independent sparse calculation uses

```text
x6 = (-1,-2,4,3,6,5), length 4
w6 = (-1,-6,3,-4,5,-2), length 15
P(x6,w6) = 1 + 6q + 11q^2 + 6q^3 + q^4 + q^5.
```

The rank-11 interval has 1,676 elements, with rank vector
`(1,12,55,140,248,339,360,287,162,59,12,1)`.
The principal lower ideal of w6 has 3,184 elements.
Every right-descent choice at the root and the independent identity
`q^d P(x,w)(q^-1) = sum_{x<=z<=w} R(x,z) P(z,w)` are checked.

The plain JSON dependency certificate contains 24,245 polynomial records:
17,267 comparable pairs (17,148 strict) and 6,978 incomparable zeros.
For every strict record it lists the descent used and every nonzero
correction term. The verifier regenerates the full correction index set,
including the omitted terms with mu=0, requires every referenced
polynomial to exist, and validates the recurrence. It independently
recomputes the R-polynomials for the root reciprocity check. It does not
compute KL polynomials by invoking either evaluator.

The even-rank n=4,...,32 bad-element table in the main JSON is only exact
arithmetic in the thesis formulas. D8 and larger bad-element KL polynomials
are not recomputed by the default command. The separate rank-4 B2 absence
test excludes a D4 FC-bottom realization of its longest principal interval;
it does not exclude models in higher D ranks.
