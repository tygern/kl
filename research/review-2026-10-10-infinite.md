# Independent referee review: infinite families

Review date: 2026-10-10. Source SHA-256: `3a39fef2c91b510f9f0ab745206efdd3afa3e2a9836a22bde28315599d8b8937`. Scope: `results/exceptional-leading.tex`, especially lines 592–948, the maximum-descent prerequisites at 265–368, and Appendix A at 1073–1078, 1126–1199. No manuscript edits were made during the review. I did not read previous AI review notes. I independently derived the main formulas, inspected the relevant supplement implementations, consulted the primary papers below, and wrote a separate arbitrary-precision scratch checker using no repository libraries.

**Subsequent revision, v0.6.0.** The manuscript now explicitly assumes a fully commutative cover in the relevant proof step and uses the five braid witnesses with length-additive translation propagation in the affine cover argument. These implement the wording repair and optional simplification recorded below. The original conclusions and line references are retained as a historical review. The independent checker was ported to standard-library Go, with `math/big` for large root samples, and checked against its original output; the preserved command is linked below.

## Verdict

I found **no substantive mathematical error or unresolved proof gap in the infinite-family results**. The reflection-cover argument works, the uniform seed is demonstrably a real root rather than merely a norm-two vector, the quadratic orbit formulas and signs are correct, the E10 lowering witness is correct, and the affine inversion count is correct. The affine finite computational input is more extensive than necessary; a short alternative below eliminates the FC catalogue from that proof. This is an optional substantial simplification, not a condition for correctness.

My assessment of correctness is high confidence for the symbolic arguments. The independently checked finite affine data and printed root witnesses agree exactly with the manuscript. I did not regenerate all 44,199 FC elements in the independent checker; the catalogue criterion and the two Go implementations were inspected, the published FC polynomial supplies the same maximum, and the 21 required base-cover exclusions were independently verified by explicit braid words without relying on the catalogue.

## Checks of the main arguments

### Reflection lemma and covers: lines 600–665

The proof of the reflection length/depth formula is valid in the stated simply laced setting. In lines 621–625, the positive integer `c` is at least one, so `c beta + (c^2-1) alpha_s` is positive. The induction in lines 628–635 gives both opposing inequalities for `ell(r_beta)=2 dp(beta)-1`. The central three-letter factor of a minimal palindromic word must be a noncommuting braid: if its outer generator equals or commutes with the middle generator, the word shortens. Thus a nonsimple reflection cannot be FC.

For the cover lemma, assume an FC cover `x`. If all common descents survive, the maximum-descent lemma forces `x=i(I)`; this has odd length, whereas a cover of a reflection has even length. If a left descent `s` is missing, the lifting property forces `x=sb`. The root `s beta` is positive, and at least two of its nonzero coordinates are off `s`, so it is not simple. Hence `r_(s beta)` is non-FC, and `bs=s r_(s beta)` is length-additive by the negative pairing. A reduced braid word for the latter reflection persists in `bs`. Taking inverses handles `sb`. Non-covers vanish by the parity proposition.

This is a useful structural argument. Its three-coordinate condition rules out the rank-two reflection `sts`, whose length-two covers are FC. The argument does not actually require terminality; the title “Covers of terminal reflections” is narrower than the stated hypotheses/proof.

One **minor logical wording repair**: line 650 should begin “Suppose that `x` is a fully commutative cover.” As written, the “Otherwise” at line 653 does not follow literally from negating the preceding conjunction: that conjunction can fail merely because `x` is non-FC while all descent inclusions hold. The intended restriction to FC covers is clear and the proof above establishes it; this is not a theorem-level gap.

### Uniform real-root seed: lines 693–725

The D_(4r) realization is correct, including the apparently delicate indefinite direction. The simple root `alpha_0=z-(1/2) sum e_i` has norm `2-r+r=2`, pairs to `-1` with `alpha_1=e_1+e_2`, and pairs to zero with the remaining D simple roots. The coefficients of the seed in the e-basis are exactly those printed. Applying `s_0` changes the z coefficient from `r-1` to one and gives e-coefficients `-1/2` at odd positions, `+1/2` at even positions. This is exactly `F alpha_0`.

Each factor `r_(e_a-e_b) r_(e_a+e_b)` changes the signs of coordinates `a,b`. The displayed product flips the `2r` even coordinates, so it belongs to the even signed-permutation group W(D_(4r)). Thus `beta_r=s_0 F alpha_0` proves real-root membership uniformly. The determinant observation is true but not needed: the nonzero orthogonal z direction and the e basis already prove linear independence for r>=3.

### Root action and terminality: lines 727–807

Direct composition gives

`N rho = (rho,delta) gamma - ((rho,gamma)+(rho,delta)) delta`,

where `N=T-1`, `N delta=0`, `N gamma=-2 delta`. Consequently `N^2 rho=-2(rho,delta)delta` and `N^3=0`, even though delta is not radical in the ambient indefinite group. The resulting binomial formula is precisely the one printed. This correctly handles the embedded affine subsystem rather than mistakenly treating its null root as globally radical.

The pairings in lines 762–765 follow directly from the tridiagonal Cartan rule, with the two exceptional positions being the branch and the node immediately beyond the embedded affine subdiagram. The special edge sums are correct, including `-ak(k-1)-rk` at edge 7–8, which is nonpositive for integer k>=0. The nonnegative first and second differences establish positivity of every coordinate.

The matching at 794–796 has exactly 2r edges, giving the upper bound 2r+1 for the independence number. The listed independent set attains it and has odd cardinality. The full-support argument using a changed row in column zero is valid. The distinctness coordinate expands to `r-1+4k+(4r-8)k^2`, as printed.

### E10: lines 908–945 and 1191–1196

All coordinates and pairings agree with the printed formula. The lowering string reduces the stated seed to alpha_9 through positive roots with strictly decreasing heights at every step. The embedded gamma and gamma+delta witnesses also work. The two relevant pairings are -3 and -1. Edge sums are nonpositive; for example the branch edge 2–9 gives zero and edge 7–8 gives `-k^2-2k`. Column one and coordinate zero prove full support and distinctness. The maximum independent set has size five. No gap found.

### Affine family: lines 822–905

I independently checked that the 27-letter word is reduced, its matrix squares to identity, its descent set is the stated I', all required adjacent column sums are positive, and its conjugation slope is exactly `d=(2,-1,1,-1,1,-1,1,-1,-1)`. The matrix formulas and support arguments are sound.

The inversion formula counts real affine roots correctly: for finite root eta, positivity starts at `m=epsilon(eta)`, and negativity of its image ends immediately below `epsilon(eta')-h-k sum eta_j d_j`. This gives precisely `max(0,A+kB)` with no missing endpoint term. Independently enumerating the 240 finite E8 roots produced these classes:

| (A,B) | Multiplicity |
|---|---:|
| (-1,-2) | 6 |
| (-1,-1) | 20 |
| (-1,0) | 1 |
| (0,-2) | 8 |
| (0,-1) | 44 |
| (0,0) | 82 |
| (0,1) | 44 |
| (0,2) | 8 |
| (1,0) | 1 |
| (1,1) | 20 |
| (1,2) | 6 |

Thus the 161 omitted contributions do vanish for all k>=0, and the length is 27+92k. Single-letter deletions give exactly 21 distinct length-26 covers, and I found and checked an explicit reduced braid word for each without consulting the FC catalogue.

## Optional simplification: remove the affine FC catalogue from the proof

This is available after the existing inversion calculation and needs only five short witnesses.

Let `C=c_0`, treat delta as a column and d as a row, and define `R=1+delta d`. Since `d delta=0`, `R^k=1+k delta d`. Since `C delta=delta`, the existing formula gives

`c_k=C R^k`.

Moreover, `R=C c_1` belongs to W, because C is an involution. The finite part of R is the identity, so in the same affine inversion formula its A-value is zero for every finite root and its B-values are exactly the previously computed B-values. Therefore

`ell(R^k)=k sum_eta max(0,B_eta)=92k`.

Together with `ell(c_k)=27+92k`, this proves that `C R^k` is length-additive. For each s in I',

`sc_k=(sC)R^k`, and `ell(sc_k)=26+92k=ell(sC)+ell(R^k)`.

The following are reduced words for the five elements sC. Brackets only identify a noncommuting braid; delete the brackets to obtain the word. Each has length 26, and the independent checker verifies both its product and reducedness.

| s | Reduced word for sC |
|---|---|
| 1 | `3[282]7564534123012856745231` |
| 3 | `1[282]7564534123012856745231` |
| 5 | `31[282]764534123012856745231` |
| 7 | `31[282]564534123012856745231` |
| 8 | `3[121]8756453423012856745231` |

Each braid persists in a reduced word for sc_k. Involutivity handles c_ks. The same maximum-descent/lifting dichotomy already printed at 891–895 then excludes every FC cover for every k>=0. This removes the need to establish the maximum FC length, to regenerate the 44,199-element catalogue, or to exclude all 21 base covers. Those computations could remain optional corroboration.

The repository's affine library already computes the translation R and verifies its length; the simplification is therefore closely aligned with existing certified data. The new ingredient presented here is the explicit five-word braid table and its propagation argument.

## Significance and scope

The strongest infinite-family contribution is the explicit construction of full-support, non-FC terminals in fixed affine/indefinite groups and uniformly in arbitrarily large rank, together with the clean reflection-cover mechanism. The seed proof and ambient unipotent calculation are particularly worthwhile structural components.

Every fixed E_n has only finitely many FC elements. This means that sufficiently long members of any family have no FC covers for the elementary length reason. **It does not by itself imply vanishing of arbitrary non-cover mu coefficients.** In this manuscript those non-covers vanish because the common descent set is maximum and has the correct parity. Thus for each fixed rank only finitely many parameter values require nontrivial cover exclusion, but the uniform lemma handles them without rank-dependent catalogues. The r=3 comparison at lines 811–814 correctly shows only that the crude length test is insufficient at the seed; it does not suggest that an FC cover exists. This is an appropriate motivation, not evidence for a counterexample.

The section does not claim a classification of terminals in infinite type or settle a full 0–1 statement for affine E8 or indefinite E_n. Its infinite vanishing statements should be evaluated with that scope in mind.

## Primary-source checks

- [Green, *On the Markov trace for Temperley–Lieb algebras of type E_n*, arXiv:0704.0283v1](https://arxiv.org/pdf/0704.0283v1): Section 3 works throughout n>=6; Proposition 4.4(iv),(v), pp. 10–11, gives precisely the cell statements used in the maximum-descent lemma. They are not restricted to finite E6–E8.
- [Biagioli–Jouhet–Nadeau, *Fully commutative elements in finite and affine Coxeter groups*](https://arxiv.org/pdf/1402.2166): Section 5.2, p. 27, displays the affine-E8 length polynomial with degree 44. The preceding discussion also records FC finiteness of every generalized E_n. Its polynomial coefficients sum to the stated count.
- [Green, *Star reducible Coxeter groups*](https://arxiv.org/abs/math/0509363): the classification includes the entire generalized E_n series. This agrees with the manuscript's explanation that noncommuting terminals in these groups cannot be FC.

## Reproduction

Run the standalone [Go checker](../tools/cmd/review-infinite/main.go) from the repository root:

```sh
go run -C tools ./cmd/review-infinite
```

The command uses only the standard library and no supplement code or input files. It checks every printed finite root witness; affine reduction, involution, descents, edge signs, inversion classes, all 21 base covers, and the five displayed braid witnesses; the translation length; and exact-integer samples of the uniform and E10 formulas through r=100 and k=1,000,000. Small fixed matrices use bounded integer arithmetic; the large root samples and their intermediate norm products use `math/big`. The parameter proofs were checked algebraically above; finite samples are cross-checks and are not substituted for those proofs.
