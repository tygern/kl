# Independent finite-theorem review — 10 October 2026

Reviewed `results/exceptional-leading.tex`, especially lines 194–545, and the finite terminal and residual-polynomial implementations. Reviewed source SHA-256: `3a39fef2c91b510f9f0ab745206efdd3afa3e2a9836a22bde28315599d8b8937`. This review did not use prior AI review conclusions. The manuscript was not edited during the review. Line references below refer to the source as read on 10 October 2026.

**Subsequent revision, v0.6.0.** The release contains the total candidate count, descriptions of both exceptional terminals as products of orthogonal reflections, and a finite-only supplement runner. These were the three optional suggestions recorded below. The original findings and line references are retained as a historical review. The independent scratch checker was ported to standard-library Go and checked against its original output; the preserved command is linked below.

## Verdict

**I found no substantive mathematical error or missing case in the finite theorem.** The argument is complete once the terminal enumeration is accepted as a computer-assisted proof; the enumeration has a valid completeness argument and its principal implementation follows that argument. The unusually important lower-endpoint lemma is a valid application of the exact results cited, not an unsupported transfer from ordinary Kazhdan–Lusztig cells. I independently verified the printed finite data and reran the substantive finite searches and residual polynomial computation described below.

Confidence in this assessment is high. It is based on mathematical inspection, source inspection, fresh executions, and a separately written direct integer-matrix check of the printed words. This is not a formal verification of Go or the compiler, nor a complete historical-priority audit.

There are **no required corrections to the finite proof identified by this review**. The optional suggestions near the end are presentation and auditability improvements, not gaps.

## Mathematical dependency audit

### Terminal reduction: lines 197–242

The descent criterion is used only after the cover case has been removed. This avoids the standard and potentially serious mistake of imposing descent inclusion on all nonzero coefficients. The induction in lines 227–242 is sound:

1. For a remaining odd gap at least three, both descent inclusions may be assumed.
2. Fully commutative elements have commuting descent sets, so the upper endpoint also has commuting descent sets.
3. If a reduced expression begins `st`, the upper endpoint is in the long position of a left rank-two string. The argument that neither `s` nor `t` is a left descent of the remaining factor `v` is correct: `tv` is reduced, while a prefix `s` of `v` would give a braid and put `t` in the upper descent set.
4. Descent inclusion and full commutativity put the lower endpoint in the same star domain; star transport preserves its full commutativity.
5. The upper length decreases by one. The lower length changes by one in either direction, so the new gap is the old gap or that gap minus two, and remains positive. Thus an incomparable pair contributes zero, a comparable pair has the correct orientation for induction, and a gap-three pair may legitimately end at a cover.

The example in lines 245–262 correctly illustrates the increasing lower star and the cover exception. No positivity theorem is secretly used: the induction transports an exact coefficient to a terminal value, a cover, or zero.

The symmetric version of star transport stated in lines 211–217 is the safe formulation. Green's original reduction and elementary identities can be checked in [Green, arXiv:0801.1650v1, §§2–3](https://arxiv.org/html/0801.1650v1). The manuscript gives a sufficient argument of its own for the finite reduction.

### Maximum common commuting descents: lines 265–317

This is correct, including the scope over generalized `E_n` and the use of a maximum independent set on the support rather than on the whole ambient graph. I checked the primary source [Green, arXiv:0704.0283v1, Proposition 4.4(iv),(v) and the end of the proof of Lemma 5.2](https://arxiv.org/html/0704.0283v1). Proposition 4.4 has exactly the right-cell characterization and singleton left/right-cell intersection claimed in the manuscript. The paper's setup is all `E_n`, `n ≥ 6`.

The proof does not need the support subgraph itself to have type `E`: `x` and `i(I)` are considered in the ambient Temperley–Lieb quotient, and the support bound only supplies `a_TL(x) ≤ |I|`. Since commuting descents provide a reduced prefix and suffix, `a_TL(x) ≥ |I|` as well. This places the two elements in the same left and right monomial cells. The singleton intersection applies and forces equality.

The distinction from ordinary KL cells in lines 272–277 is essential and correct. The maximum-versus-maximal warning in lines 319–323 is also correct. There is no circular use of the desired KL bound in this lemma.

### Parity and the finite support data: lines 325–443

For an eligible non-cover, the subword property gives the required support containment and the maximum-descent lemma forces the commuting lower endpoint. The parity conclusion then follows directly from the definition of `mu`. Covers have already been handled separately.

All stated matching bounds in lines 433–439 are valid: a matching gives an upper bound `|support| − |matching|` for an independent set, whether or not one invokes the stronger bipartite matching theorem. The displayed independent descent set attains each bound. In particular, the length-50 `E8` terminal really has support independence number four; no assumption about all maximal independent sets having the same size is involved.

### Terminal geometry and parabolic pruning: lines 355–368 and 450–499

For a right descent `s`, a reduced word ending with the noncommuting pair `ts` exists exactly when `t` is a right descent of `ws`. The identity `(ws)(alpha_t)=w(alpha_s+alpha_t)` therefore proves the printed terminal test. Testing the inverse gives precisely the left condition.

The pruning lemma is sufficient and correctly one-sided. It does not erroneously retain only two-sided terminals at intermediate stages. If `w=av` is length additive and right-terminal, every reduced suffix word for `v` can be appended to a reduced word for `a`, so `v` must be right-terminal in the smaller parabolic.

The coset queue is exhaustive: stripping right `J`-descents finds the unique minimal representative of the same right coset; left multiplication by generators acts on the right-coset space and reaches every coset from the identity. The induction invariant in lines 491–495 then proves that every right-terminal element is tested. The last inverse test imposes the other side. Group orders are useful independent checks, but neither the proof nor the code relies only on a declared `complete: true` flag or an expected count to supply exhaustiveness.

The recursion is therefore a genuine exhaustive search with a proved pruning rule, rather than a search up to an unexplained length bound. The internal assertion limiting a reduced E8 word to length 120 is consistent with the finite root system's 120 positive roots; it is a failure guard, not a truncation that silently discards candidates.

### Final finite theorem and products: lines 512–545

The commuting-terminal case is correct by the subword property and the descent criterion. The only remaining odd-gap table entries are the two ambient copies of the same `D6` pair. Standard parabolic invariance of ordinary KL polynomials applies to that pair. [Gern, arXiv:1304.6074, Lemma 4.5.5](https://arxiv.org/pdf/1304.6074) explicitly gives `mu(x6,w6)=1`; I also reran its independent ambient calculation.

The direct-product argument is correct. For an odd total gap with at least two nontrivial components, summing the component degree bounds gives a strictly smaller degree than the degree defining `mu`. Thus products cannot create a new value above one. The proof is not implicitly multiplying leading coefficients at incompatible degrees.

## Implementation audit

The core code read was:

- `supplement/internal/parabolic/parabolic.go`, especially root closure and packing (101–185), right/left multiplication (197–214), descent and terminal tests (217–258), minimal cosets (360–390), composition (408–434), and diagram orders/chains (591–715).
- `supplement/cmd/terminals-recursive/main.go`, especially the extension and chain loops (127–212) and final inverse test/classification (224–239).
- `supplement/cmd/e8-d7/main.go`, including full `D7` BFS (279–307), the fundamental-weight orbit coset construction (309–362), and the complete product/filter loop (410–430).
- `supplement/cmd/matrix-search/e7.go`, whose terminal test uses actual descent transitions after multiplication, providing a useful alternative to the root-sum implementation.
- `supplement/cmd/d6-ambient/main.go`, particularly full subword-ideal construction and the KL recurrence (354–535).

The packed representation retains the images of all simple roots, including those outside the current parabolic, so distinct ambient elements cannot collapse through an insufficient subsystem encoding. E8's 240 roots fit in the eight-bit root identifiers. Products are composed in the correct order, and inverse words are correctly reversed.

In the separate D7 program, the orbit vector has pairing one with node 0 and zero with every D7 node. Its orbit supplies the expected 2160 distinct cosets, and the code explicitly checks no right D7 descents for the representatives. The complete `2160 × 280` candidate product set is checked. This is more independent than simply repeating the same recursive chain.

The residual-polynomial code uses the correct recurrence, including the exponent on each correction term, and computes all polynomials needed inside the complete lower ideal. It does not read the desired polynomial and merely echo it: the saved D6 certificate is only an additional cross-check after the computation.

## Fresh checks performed

All executions returned success. Temporary JSON outputs are under `/tmp/kl-finite-review-*`.

| Check | Result |
|---|---|
| Full E6 root-index enumeration | Order 51,840; 22 commuting terminals; one noncommuting terminal of length 7 |
| Full E7 integer-matrix enumeration | Order 2,903,040; 36 commuting terminals; noncommuting lengths 7, 8, 15, 28 |
| Recursive E7 classification | 576 right-terminals; same four noncommuting rows, descent masks and eligible-bottom data |
| Four recursive E8 chains | Every chain gives the same 2160 right-terminals and 64 terminals; 58 commuting and 6 noncommuting |
| Separate E8/D7 enumeration | Same six noncommuting terminals, lengths 7, 8, 8, 15, 28, 50 |
| D6 ambient KL computation | `1+6q+11q²+6q³+q⁴+q⁵`, `mu=1`; lower ideal 3184 and interval 1676 in E7, E8 and D6 |
| Relevant Go tests | `internal/parabolic`, `cmd/terminals-recursive`, `cmd/e8-d7`, `cmd/d6-ambient` all pass |

I additionally wrote and executed a small direct integer-column calculation separate from the supplement. Its self-contained checks are now preserved as [tools/cmd/review-finite/main.go](../tools/cmd/review-finite/main.go), a standalone Go command using only the standard library. It generates the finite roots by reflection closure, computes word length by counting inverted positive roots, and checks all eleven table words. Every word has its printed length, equals its inverse, passes terminality on both sides, and has exactly the displayed common descent set. Brute enumeration of independent subsets of each support confirms maximal cardinality. It also checks the orthogonal-reflection descriptions below. During the original review, the scratch checker's matrices matched all noncommuting E8 words from every fresh chain and the separate D7 enumeration, and matched the fresh E6/E7 classifications. This directly addresses transcription, representation convention and reducedness risks. Those cross-search comparisons are a historical record; the preserved Go command checks the table and descriptors without reading certificates or depending on temporary files. It does not replace the exhaustive enumeration commands.

Representative reproduction commands, from the repository root:

```sh
go run -C tools ./cmd/review-finite
GOCACHE=/tmp/kl-finite-review-gocache go run -C supplement ./cmd/terminals-recursive -rank 8 -chains-certificate /tmp/kl-finite-review-e8-chains.json
GOCACHE=/tmp/kl-finite-review-gocache go run -C supplement ./cmd/e8-d7
GOCACHE=/tmp/kl-finite-review-gocache go run -C supplement ./cmd/matrix-search -rank 7
GOCACHE=/tmp/kl-finite-review-gocache go run -C supplement ./cmd/enumerate-bad -rank 6
GOCACHE=/tmp/kl-finite-review-gocache go run -C supplement ./cmd/terminals-recursive -rank 7
GOCACHE=/tmp/kl-finite-review-gocache go run -C supplement ./cmd/d6-ambient -out /tmp/kl-finite-review-d6.json -d6cert ../results/d6-recurrence-certificate.json
```

## Significance and optional improvements

The finite result is a meaningful completion of the finite simply laced case, assuming the priority claim survives the separate literature audit. Its substance is the compact terminal classification: two exceptional obstructions beyond the embedded type-D list, both eliminated without new KL computation. This is stronger and more explanatory than a raw table of all coefficients. The parabolic pruning turns the enormous E8 group into only 138,240 last-stage candidates on the main chain, with 143,660 candidates over all stages, and the manuscript explains why that reduction is exhaustive.

The conceptual machinery is largely inherited: descent/star transport, the monomial-cell singleton argument, parabolic factorization and the D6 coefficient. The manuscript already gives those attributions accurately in lines 128–141. I would describe the new finite contribution as a rigorous computational classification with an efficient structural deduction, rather than as a fundamentally new general theory of KL coefficients. That is a scope assessment, not a criticism of correctness or publishability.

Three optional improvements would make referee work easier:

1. **Auditability, low severity:** near lines 501–507, print the total main-chain candidate count 143,660, or a compact stage-count table. The source and certificates contain it, but the number makes the actual scale of the proof especially clear.
2. **Reproduction ergonomics, low severity:** add a finite-only runner target. The all-results runner described at lines 1211–1228 also verifies extensive infinite-family and supplementary calculations. A referee assessing Theorem 1.1 could benefit from one command for just terminal exhaustiveness, table matrices, and the residual D6 value. Existing individual commands already suffice, so this is not a missing verification.
3. **Conceptual presentation, low severity:** bring one compact intrinsic description of each new terminal into the paper, rather than leaving every such descriptor in the certificate. Lines 565–568 already announce orthogonal-reflection decompositions. In the paper's coordinates, the certificate gives the E8 terminal as `r_(1,2,3,3,2,2,1,2) r_(1,3,4,3,2,1,0,2)`. The E7 terminal uses the four orthogonal roots `(0,0,1,1,1,1,1)`, `(0,1,1,1,1,0,1)`, `(0,1,2,1,0,0,1)`, `(1,2,2,2,1,1,1)`. I separately checked pairwise orthogonality, norm two, and equality of the resulting reflection products to the printed words. A short display of this kind would make the two new objects easier to recognize and remember. It would supplement, not replace, the exhaustive classification proof.

None of these suggestions is needed to establish the finite theorem. On the finite part alone, I would support publication in a specialist algebraic-combinatorics or representation-theory venue after ordinary editorial revision and preservation of an accessible, fixed computational supplement. I would not recommend rejection or a demand for a wholly noncomputational proof: the mathematical compression and verified exhaustive computation are legitimate contributions. I also would not sell the finite proof as a new general mechanism beyond the existing theory. I do not recommend manufacturing a mathematical objection where this review found none.
