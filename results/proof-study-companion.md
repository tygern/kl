# Studying the proofs of the exceptional leading coefficient theorem

The proof of the [exceptional leading coefficient theorem](exceptional-leading.tex) uses star reduction, maximum common descents, and exhaustive parabolic pruning. The infinite families use the same descent argument, explicit root constructions, and a separate proof excluding fully commutative covers. The notes below include examples and exercises for these arguments.

Work through the exercises before reading the solution sketches, identifying where each hypothesis is used. References use result names and source labels, which are independent of theorem numbering.

## 1 The finite proof and its dependencies

For an odd gap $d=\ell(w)-\ell(x)>0$, the coefficient in question is

$$
\mu(x,w)=[q^{(d-1)/2}]P_{x,w}(q).
$$

It is zero by definition for even gaps and for $x\not<w$. A cover has gap one and coefficient one. The phrase “leading coefficient” here means the coefficient at the largest *allowed* degree; that coefficient can vanish even when the polynomial is nonzero.

The finite proof has the following dependencies.

| Argument | Conclusion | Required results |
|---|---|---|
| Descent criterion and terminal reduction | It suffices to check terminal upper endpoints and covers | Standard KL descent identities and star transport; preservation of full commutativity |
| Computational classification | A complete finite list of noncommuting terminals | Faithful group arithmetic, the geometric terminal test, and exhaustive parabolic pruning |
| Maximum common commuting descents | Each listed terminal has one eligible FC lower endpoint | Maximum size of the descent set and the two quoted Fan–Green monomial cell facts |
| Parity and the inherited coefficient | The remaining values are zero or one | Table lengths and Gern’s $D_6$ value |

The classification requires both program verification and a proof of completeness. The maximum-descent lemma follows from the two cited results in monomial cell theory.

Compare *Maximum common commuting descents* with Gern’s Lemma 4.5.8, and *Terminal reduction* with the reduction in Gern’s Theorem 4.5.11. For the finite exceptional types, the reduction is applied to the classified terminals.

**Exercise 1.** Explain why the coefficient vanishes by parity for the length-28 terminal with lower endpoint of length four, while the length-15 terminal with lower endpoint of length four requires a coefficient calculation. Which coefficient is needed in the latter case?

## 2 Follow one star operation completely

Read *Star reduction*, especially the definition of $D_L(s,t)$ and the transport identity (`eq:star`). A star operation acts on an element with exactly one of $s,t$ as a left descent. Multiplication by exactly one of these two generators leaves the element in that domain. This partner can have either greater or smaller length.

Use the path $0-1-2$, with $s_i$ the adjacent transposition of positions $i+1$ and $i+2$ in $S_4$. Products act on the left; their rightmost factor acts first. Right multiplication swaps positions in one-line notation. Consider

$$
x=s_1,\qquad w=s_1s_0s_2s_1.
$$

For the pair of generators $(s_1,s_0)$, the relevant data are:

| Element | One-line permutation | Length | Left descents |
|---|---|---:|---|
| $x=s_1$ | $1324$ | 1 | $\{s_1\}$ |
| $w=s_1s_0s_2s_1$ | $3412$ | 4 | $\{s_1\}$ |
| ${}^*x=s_0s_1$ | $2314$ | 2 | $\{s_0\}$ |
| ${}^*w=s_0s_2s_1$ | $2413$ | 3 | $\{s_0,s_2\}$ |

All four elements have exactly one left descent in $\{s_1,s_0\}$. The star therefore increases the lower length and decreases the upper length. The new pair is comparable: deleting $s_2$ from $s_0s_2s_1$ gives $s_0s_1$. It is a cover, so

$$
\mu(s_1,s_1s_0s_2s_1)=\mu(s_0s_1,s_0s_2s_1)=1.
$$

The generator $s_2$ is a left descent of the new upper endpoint but not of the new lower endpoint. Since the lower endpoint equals $s_2{}^*w$, the pair satisfies the cover exception in the descent criterion.

The general induction uses the upper length as its decreasing quantity. If the old gap is an odd $d\geq3$, the new gap is $d$ or $d-2$, both positive. Transport is stated using the symmetric coefficient $\widetilde\mu$ because the star partners need not be comparable. If they are incomparable, the transported value is zero; if comparable, the lower star partner has smaller length, as required for induction. The transport identity is an equality of $\mu$ coefficients.

**Exercise 2.** Recover each length in the table by counting inversions. Recover each left descent set by counting descents of the inverse permutation. Check the two alternative multiplications in the definition of each star.

**Exercise 3.** In the general proof, prove that $x$ belongs to $D_L(s,t)$ when descent inclusion holds, $x$ is fully commutative, and the chosen reduced prefix of $w$ is $st$. Identify exactly where $d\geq3$ is used.

## 3 Understand the maximum-descent argument

Read *Maximum commuting descents* (`lem:maxI`). Let $K$ contain the support of the FC element $x$, and let $I$ be a maximum independent set in the induced graph on $K$. Suppose $I\subseteq L(x)\cap R(x)$.

Since the left descents in $I$ commute, $x=i(I)v$ length-additively. Since they are also right descents, $x=u\,i(I)$ length-additively. The prefix and suffix can overlap; a factorization $i(I)z\,i(I)$ with disjoint copies of $I$ need not exist.

The number

$$
a_{\rm TL}(x)=\max\{|A|:x=u\,i(A)\,v\text{ length-additively},\ A\text{ independent}\}
$$

is the maximum number of generators in a commuting factor of a reduced expression. Since $i(I)$ is a reduced prefix, $a_{\rm TL}(x)\geq|I|$. Every commuting factor uses an independent subset of $K$. Since $I$ has maximum cardinality, the reverse inequality holds. Hence $a_{\rm TL}(x)=|I|$.

By the quoted Fan–Green characterization, $x$ belongs to the monomial right cell of $i(I)$ and to its monomial left cell. These two cells intersect in one element, and $i(I)$ is already in the intersection. Thus $x=i(I)$.

The equality uses the singleton intersection theorem for cells of the Temperley–Lieb monomial basis. That theorem concerns monomial cells, which differ from ordinary KL cells. For its proof, see Green’s Proposition 4.4(iv)–(v) and Fan’s corresponding results cited in the manuscript.

On the path $0-1-2$, the independent set $\{1\}$ is maximal under inclusion but is not maximum: $\{0,2\}$ is larger. The FC element $s_1s_0s_2s_1$ has both descent sets $\{s_1\}$ and contains the commuting factor $s_0s_2$. Consequently its $a_{\rm TL}$ value is two, not one. This is why “maximum” cannot be weakened to “maximal.”

For a terminal $b$, eligibility includes $x\leq b$. By the subword property, $\mathrm{supp}(x)\subseteq\mathrm{supp}(b)$, so we can take $K$ to be the support of $b$. With $x=i(I)$, the coefficient vanishes whenever the gap $\ell(b)-|I|$ is even.

**Exercise 4.** Write the maximum-descent proof in five lines, marking the two lines that use the cited cell results. Explain why the support assumption is needed for the upper bound on $a_{\rm TL}(x)$.

## 4 Prove that the enumeration cannot miss an element

Read *Exhaustiveness* (`sec:exhaustive`) and the parabolic pruning lemma (`lem:pruning`). Fix $J\subset K$. Every $w\in W_K$ has a unique decomposition

$$
w=av,\qquad v\in W_J,\qquad R(a)\cap J=\varnothing,
\qquad \ell(w)=\ell(a)+\ell(v).
$$

If a reduced word of $v$ ends in noncommuting generators, prefixing it with a reduced word of $a$ produces a reduced word of $w$ with the same noncommuting suffix. Thus every right-terminal $w$ has a right-terminal parabolic factor $v$. The converse is false: the product must still pass the terminal test in $K$.

At stage $K$, the retained set contains exactly the right-terminal elements of $W_K$. The proof has two parts:

- Soundness: every retained candidate passes the right-terminal test.
- Completeness: a right-terminal $w$ has a unique parabolic factorization. By the pruning lemma, its factor $v$ is in the preceding retained set; its representative $a$ is in the complete coset list. Therefore the algorithm tests $w$.

Representatives are enumerated by a queue. From a known minimal representative $a$, left multiply by a generator of $K$, then strip right $J$-descents. Stripping decreases length and reaches the unique minimal representative of the same coset. Left multiplication acts transitively on the cosets, so every coset occurs when the queue is exhausted. The program also compares the number of representatives with the group index. For these finite parabolics, the procedure terminates without a search-depth cutoff.

### A complete example in type A2

Let $K=\{s,t\}$ with $m(s,t)=3$, and $J=\{s\}$. The preceding retained set is $\{e,s\}$, and the minimal representatives of $W_K/W_J$ are $\{e,t,st\}$. Their products are:

| Representative $a$ | Factor $v$ | Product | Right-terminal? |
|---|---|---|---|
| $e$ | $e$ | $e$ | Yes |
| $e$ | $s$ | $s$ | Yes |
| $t$ | $e$ | $t$ | Yes |
| $t$ | $s$ | $ts$ | No |
| $st$ | $e$ | $st$ | No |
| $st$ | $s$ | $sts$ | No |

Thus the retained set is $\{e,s,t\}$. All three pass the inverse test as well. The remaining products fail the right-terminal test, although every factor $v$ is right-terminal.

In $E_8$, the final stage combines 240 representatives with 576 right-terminals of the $E_7$ parabolic. It tests 138,240 products, retaining 2,160 right-terminals. The inverse test leaves 64 terminals: 58 commuting products and six noncommuting elements. These counts concern the final stage; the total over all stages is 143,660 tested products.

**Exercise 5.** Reproduce the representative queue for the A2 example. Explain why stripping right descents preserves the right coset, whereas choosing arbitrary left descents would not justify the same algorithm.

**Exercise 6.** Write the completeness induction in your own words without mentioning a program. Then identify where the actual implementation must provide exact group equality, exact descent tests, and complete coset enumeration.

For implementation study, begin with `MinimalCosets` in [parabolic.go](https://github.com/tygern/kl/blob/v0.5.0/supplement/internal/parabolic/parabolic.go), then `extendRightTerminals` in [terminals-recursive](https://github.com/tygern/kl/blob/v0.5.0/supplement/cmd/terminals-recursive/main.go). These functions implement the two parts of the argument above. The links refer to fixed release sources and are accessible from the proof archive. Separate enumeration programs also check the output lists.

## 5 Build a terminal reflection from its root pairings

Read *Reflections and their covers*, then *The rank-uniform construction*. For a positive real root $\beta$ of norm two, set $m_j=(\alpha_j,\beta)$. The reflection formula is

$$
r_\beta(\alpha_j)=\alpha_j-m_j\beta.
$$

Suppose $\beta$ has full support on more than one node. If $m_j>0$, then the integral pairing is at least one, and the reflected root is negative. If $m_j\leq0$, the reflected root is positive. Thus $j$ is a right descent exactly when $m_j>0$. Since a reflection is an involution, the left descent set is the same.

For an edge $s\sim t$ incident to a descent with $m_s+m_t\leq0$,

$$
r_\beta(\alpha_s+\alpha_t)
=\alpha_s+\alpha_t-(m_s+m_t)\beta>0.
$$

This is the geometric terminal test. We therefore construct real roots whose pairings satisfy these sign and edge conditions. In indefinite type, a vector of norm two need not be a real root. The seed belongs to a Weyl-group orbit of a simple root, and each iterate is obtained by applying reflections in real roots.

### Why the orbit formula is quadratic

Put $T=r_\gamma r_{\gamma+\delta}$ with $(\gamma,\gamma)=2$ and $(\gamma,\delta)=(\delta,\delta)=0$. If $N=T-1$, composing the reflections gives

$$
N\rho=(\rho,\delta)\gamma-
\big((\rho,\gamma)+(\rho,\delta)\big)\delta.
$$

Consequently $N\delta=0$, $N\gamma=-2\delta$, and $N^3=0$. For $(\beta,\delta)=-a$ and $(\beta,\gamma)=-c$,

$$
N\beta=-a\gamma+(a+c)\delta,\qquad N^2\beta=2a\delta.
$$

The binomial expansion terminates:

$$
T^k\beta=\beta+kN\beta+\binom{k}{2}N^2\beta
=\beta-ak\gamma+(ak^2+ck)\delta.
$$

The null vector belongs to an embedded affine subsystem. It need not be orthogonal to the entire ambient root space; the coefficient of $k^2\delta$ is $a=-(\beta,\delta)$.

### One further E13 calculation

For $r=3$, the seed is

$$
\beta_{3,0}=(2,3,5,5,4,4,3,3,2,2,1,1,3).
$$

Here $a=2$ and $c=3$, so $\beta_{3,1}=\beta_{3,0}-2\gamma+5\delta$. Explicitly,

$$
\beta_{3,1}=(10,19,29,26,20,17,11,8,2,2,1,1,16),
$$

with simple-root pairings

$$
(m_0,\ldots,m_{12})=(1,-1,-3,3,-3,3,-3,3,-6,1,-1,1,3).
$$

The positive entries still occur at $I_3=\{0,3,5,7,9,11,12\}$. The exceptional edge sums at $7-8$ and $8-9$ are now $-3$ and $-5$. The required inequalities hold for every nonnegative integer $k$ by the general pairing formulas.

### Why vanishing needs another argument for covers

A reflection has odd length. If its maximum common descent set has odd size, the only eligible FC non-cover has an even gap, so its coefficient vanishes. Covers have coefficient one, so they require a separate exclusion argument.

Suppose $x$ is an FC cover of $b$. If all descents of $b$ are descents of $x$ on both sides, then $x=i(I)$ by the maximum-descent lemma, contradicting the parity of a cover’s length. Otherwise, by the lifting property, $x=sb$ or $x=bs$ for a descent $s$ of $b$. Each of these elements has a reduced word containing a braid, as proved in the reflection argument, contradicting full commutativity.

When $\beta$ has at least three nonzero coordinates, $s\beta$ is still nonsimple: only the coordinate at $s$ changes, so at least two others remain nonzero.

The affine family has a different upper endpoint, an involution that is not a reflection. Its cover exclusion uses a reduced-prefix argument. Write $C=c_0$ and $R=1+\delta'd$, with $\delta'$ a column and $d$ a row. By the column formula, $c_k=CR^k$. Counting inversions, we obtain $\ell(R^k)=92k$ and $\ell(c_k)=27+92k$. Thus the product is length-additive. For each $s\in I'$, the factorization $sc_k=(sc_0)R^k$ is also length-additive, with length $26+92k$.

The five displayed reduced words for $sc_0$ each contain a noncommuting braid. After appending a reduced word for $R^k$, the resulting reduced expression for $sc_k$ contains the same braid. The inverse $c_ks$ is also non-FC. Thus no FC cover exists, by the lifting argument above. Separate programs check the FC catalogue and all 21 base covers as supplementary verification.

For each fixed $E_n$, the FC elements form a finite set. Elements of sufficiently large length therefore have no FC covers. For the constructed families, no FC covers exist at any parameter by the reflection and translation arguments. Coefficients at non-covers vanish by maximum common descents and parity.

**Exercise 7.** Derive the formula for $N\rho$ directly from the two reflections, keeping their order. Then derive $N^3=0$ and the orbit formula without assuming that $\delta$ is in the ambient radical.

**Exercise 8.** Recompute $m_2$ and $m_{12}$ from $\beta_{3,1}$ using the branch at node 2. Use the matching printed in the theorem proof to show that $I_3$ is maximum. Explain why the coefficients at non-covers vanish and why a separate cover lemma is needed.

## 6 Solution sketches

1. The gaps are $28-4=24$ and $15-4=11$. The first is even, so $\mu=0$ by definition. The second requires $[q^5]P_{x_6,w_6}$, which is one in the displayed $D_6$ polynomial.

2. The inversion counts are $1,4,2,3$. The inverses are respectively $1324,3412,3124,3142$; their descents give the table’s left descents. For $x$, multiplication on the left by $s_1$ gives $e$, outside the star domain, whereas $s_0x=s_0s_1$ is inside. For $w$, $s_1w=s_0s_2s_1$ is inside; $s_0w=s_0s_1s_0s_2s_1=s_1s_0s_1s_2s_1$ has both $s_0$ and $s_1$ as left descents, so is outside.

3. The initial $s$ is a descent of $w$ and hence of $x$. An FC element cannot have noncommuting $s,t$ as simultaneous descents, so $t$ is not a descent of $x$. When $d\geq3$, both possible new gaps, $d$ and $d-2$, are positive. The gap-one case must already have been handled.

4. Since $i(I)$ is a commuting prefix, the lower bound holds. Since every commuting factor has support in $K$ and $I$ has maximum cardinality, the upper bound holds. By the right-cell characterization and its left analogue, $x$ belongs to the two cells of $i(I)$. Their intersection is the singleton containing $i(I)$, so $x=i(I)$. Without support containment, a larger commuting factor could use generators outside $K$, beyond the scope of the maximum-cardinality bound.

5. From $e$, left multiplication by $s$ strips back to $e$, and multiplication by $t$ yields $t$. From $t$, multiplication by $s$ yields $st$; multiplication by $t$ yields $e$. From $st$, multiplication by $s$ yields $t$; multiplication by $t$ gives $tst=sts$, which strips on the right to $st$. Right multiplication by an element of $W_J$ preserves the coset $aW_J$. Left multiplication generally changes it.

6. Start from $\{e\}$. At each stage factor any right-terminal $w$ as $av$. By the pruning lemma and the induction hypothesis, $v$ belongs to the preceding list. The complete coset list contains $a$, and the exact terminal test accepts $w$. Conversely every accepted product satisfies the definition. Exact equality tests are needed to detect repeated representatives and products, and exact descents are needed both for stripping and for the terminal test.

7. Set $u=(\rho,\gamma)$ and $v=(\rho,\delta)$. The first reflection sends $\rho$ to $\rho-(u+v)(\gamma+\delta)$, whose pairing with $\gamma$ is $-u-2v$. Reflecting in $\gamma$ therefore gives $\rho+v\gamma-(u+v)\delta$. The image of $N$ lies in the span of $\gamma,\delta$; $N$ maps that span into the line through $\delta$ and maps that line to zero. Substitution gives the displayed binomial formula.

8. At the branch, $m_2=2(29)-19-26-16=-3$ and $m_{12}=2(16)-29=3$. The matching $(0,1),(2,12),(3,4),(5,6),(7,8),(9,10)$ has six edges, so an independent set uses at most one endpoint of each edge and the one unmatched node: at most seven nodes. The displayed $I_3$ attains seven. Both the reflection length and seven are odd, so the eligible gap is even. A hypothetical FC cover would have coefficient one, which is why the separate cover exclusion remains essential.

## 7 Review questions

Explain why the coefficient is invariant under star reduction, why an eligible FC endpoint is unique under the maximum-descent hypothesis, why every terminal occurs in the parabolic enumeration, and which coefficient is computed in Gern’s thesis.

For the infinite families, explain why the seed is a real root, why the terminal inequalities hold for every parameter, and why no FC cover exists. Identify each finite computational input and each symbolic argument. The dependencies are listed in Appendix A.
