# Studying the proofs of the exceptional leading coefficient theorem

This companion is for working through [the manuscript](exceptional-leading.tex) with a background in mathematics and software engineering. It develops four arguments: star reduction, maximum common descents, exhaustive parabolic pruning, and the root construction. The finite theorem is the first destination; the infinite families reuse its descent argument but require a separate treatment of covers.

Read one part at a time, do its exercises before reading the solution sketches, and explain the argument aloud without the paper. The goal is to identify exactly where each hypothesis enters. References below use named results and source labels so they remain usable if theorem numbering changes.

## 1 The finite proof and its dependencies

For an odd gap $d=\ell(w)-\ell(x)>0$, the coefficient in question is

$$
\mu(x,w)=[q^{(d-1)/2}]P_{x,w}(q).
$$

It is zero by definition for even gaps and for $x\not<w$. A cover has gap one and coefficient one. The phrase “leading coefficient” here means the coefficient at the largest *allowed* degree; that coefficient can vanish even when the polynomial is nonzero.

The finite proof has the following dependencies.

| Argument | What it establishes | What it uses |
|---|---|---|
| Descent criterion and terminal reduction | It suffices to handle terminal upper endpoints, with covers handled separately | Standard KL descent identities and star transport; preservation of full commutativity |
| Computational classification | A complete finite list of noncommuting terminals | Faithful group arithmetic, the geometric terminal test, and exhaustive parabolic pruning |
| Maximum common commuting descents | Each listed terminal has one eligible FC lower endpoint | Maximum size of the descent set and the two quoted Fan–Green monomial cell facts |
| Parity and the inherited coefficient | The remaining values are zero or one | Table lengths and Gern’s $D_6$ value |

The classification is a computational theorem with a mathematical completeness argument. The maximum-descent lemma is an abstract consequence of imported cell theory. Keeping these two dependencies distinct helps you decide what to check in a program and what to study in a source proof.

For your thesis connection, first compare the manuscript’s *Maximum common commuting descents* with Gern’s Lemma 4.5.8, and *Terminal reduction* with the reduction used in Gern’s Theorem 4.5.11. The finite exceptional extension uses this familiar reduction strategy; the new classification provides the cases to which it applies.

**Exercise 1.** Explain why parity settles the length-28 terminal with lower endpoint of length four, but does not settle the length-15 terminal with lower endpoint of length four. Which coefficient is needed in the latter case?

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

Notice that $s_2$ is a left descent of the new upper endpoint but not of the new lower endpoint. This does not force zero: the lower endpoint equals $s_2{}^*w$, precisely the cover exception in the descent criterion.

The general induction uses the upper length as its decreasing quantity. If the old gap is an odd $d\geq3$, the new gap is $d$ or $d-2$, both positive. Transport is stated using the symmetric coefficient $\widetilde\mu$ because the star partners need not be comparable. If they are incomparable, the transported value is zero; if comparable, the positive gap puts them in the correct order for induction. The identity transports $\mu$; no equality of entire KL polynomials is being asserted.

**Exercise 2.** Recover each length in the table by counting inversions. Recover each left descent set by counting descents of the inverse permutation. Check the two alternative multiplications in the definition of each star.

**Exercise 3.** In the general proof, explain why descent inclusion and full commutativity put $x$ in $D_L(s,t)$ whenever the chosen reduced prefix of $w$ is $st$. Identify exactly where $d\geq3$ is used.

## 3 Understand the maximum-descent argument

Read *Maximum commuting descents* (`lem:maxI`). Let $K$ contain the support of the FC element $x$, and let $I$ be a maximum independent set in the induced graph on $K$. Suppose $I\subseteq L(x)\cap R(x)$.

Commuting left descents can be placed together at the start of a reduced word, giving $x=i(I)v$ length-additively. Commuting right descents similarly give $x=u\,i(I)$. These are two factorizations of the same element. They do not assert a factorization $i(I)z\,i(I)$ with two disjoint copies of $I$.

The number

$$
a_{\rm TL}(x)=\max\{|A|:x=u\,i(A)\,v\text{ length-additively},\ A\text{ independent}\}
$$

measures the largest commuting factor that can occur in a reduced expression. The prefix $i(I)$ gives $a_{\rm TL}(x)\geq|I|$. Every such factor uses generators in $K$, so maximum cardinality gives the reverse inequality. Hence $a_{\rm TL}(x)=|I|$.

Here is the point that requires cell theory. The quoted Fan–Green characterization puts $x$ in the monomial right cell of $i(I)$ and in its monomial left cell. These two cells intersect in one element, and $i(I)$ is already in the intersection. Thus $x=i(I)$.

The inequalities alone do not finish the argument. The singleton intersection is the imported rigidity statement. These are cells of the Temperley–Lieb monomial basis; ordinary KL cells cannot simply be substituted. To study this dependency further, follow the manuscript’s precise citations to Green’s Proposition 4.4(iv)–(v) and Fan’s corresponding results. This companion explains their application, not a new proof of those cell theorems.

On the path $0-1-2$, the independent set $\{1\}$ is maximal under inclusion but is not maximum: $\{0,2\}$ is larger. The FC element $s_1s_0s_2s_1$ has both descent sets $\{s_1\}$ and contains the commuting factor $s_0s_2$. Consequently its $a_{\rm TL}$ value is two, not one. This is why “maximum” cannot be weakened to “maximal.”

For a terminal $b$, eligibility includes $x\leq b$. The subword property gives $\mathrm{supp}(x)\subseteq\mathrm{supp}(b)$, which supplies the $K$ needed above. Once $x=i(I)$ is forced, the gap $\ell(b)-|I|$ determines whether parity disposes of the coefficient.

**Exercise 4.** Write the maximum-descent proof in five lines, marking the two lines that invoke imported cell results. Explain why the support assumption is needed for the upper bound on $a_{\rm TL}(x)$.

## 4 Prove that the enumeration cannot miss an element

Read *Exhaustiveness* (`sec:exhaustive`) and the parabolic pruning lemma (`lem:pruning`). Fix $J\subset K$. Every $w\in W_K$ has a unique decomposition

$$
w=av,\qquad v\in W_J,\qquad R(a)\cap J=\varnothing,
\qquad \ell(w)=\ell(a)+\ell(v).
$$

If a reduced word of $v$ ends in noncommuting generators, prefixing it with a reduced word of $a$ gives the same defect in $w$. Thus every right-terminal $w$ has a right-terminal parabolic factor $v$. The converse is false: the product must still pass the terminal test in $K$.

The invariant in the paper’s pseudocode is that the retained set at stage $K$ contains *all and only* the right-terminal elements of $W_K$. Its two directions have different jobs:

- Soundness: every retained candidate passes the right-terminal test.
- Completeness: a right-terminal $w$ has a unique parabolic factorization, and pruning puts its factor $v$ in the preceding retained set; its representative $a$ is in the complete coset list. Therefore the algorithm tests $w$.

Representatives are enumerated by a queue. From a known minimal representative $a$, left multiply by a generator of $K$, then strip right $J$-descents. Stripping decreases length and reaches the unique minimal representative of the same coset. Left multiplication connects all cosets, so closure of the queue proves exhaustiveness. The group index provides an additional count check. The procedure is applied to finite parabolics here; it has no arbitrary search-depth cutoff.

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

Thus the retained set is $\{e,s,t\}$. All three pass the inverse test as well. The table also shows why keeping every product would be wrong even though every factor $v$ is right-terminal.

In $E_8$, the final stage combines 240 representatives with 576 right-terminals of the $E_7$ parabolic. It tests 138,240 products, retaining 2,160 right-terminals. The inverse test leaves 64 terminals: 58 commuting products and six noncommuting elements. These are counts for the final stage, not the total work of constructing earlier stages.

**Exercise 5.** Reproduce the representative queue for the A2 example. Explain why stripping right descents preserves the right coset, whereas choosing arbitrary left descents would not justify the same algorithm.

**Exercise 6.** Write the completeness induction in your own words without mentioning a program. Then identify where the actual implementation must provide exact group equality, exact descent tests, and complete coset enumeration.

For implementation study, begin with `MinimalCosets` in [parabolic.go](https://github.com/tygern/kl/blob/v0.5.0/supplement/internal/parabolic/parabolic.go), then `extendRightTerminals` in [terminals-recursive](https://github.com/tygern/kl/blob/v0.5.0/supplement/cmd/terminals-recursive/main.go). These functions implement the two parts of the argument above. The links identify the release sources and also work when this companion is read from the proof archive. Other engines provide cross-checks; reading all of them is not a prerequisite for understanding completeness.

## 5 Build a terminal reflection from its root pairings

Read *Reflections and their covers*, then *The rank-uniform construction*. For a positive real root $\beta$ of norm two, set $m_j=(\alpha_j,\beta)$. The reflection formula is

$$
r_\beta(\alpha_j)=\alpha_j-m_j\beta.
$$

When $\beta$ has full support on more than one node, positivity of $m_j$ makes this root negative: the integral positive pairing is at least one, and every coordinate of $\beta$ is positive. Nonpositive $m_j$ makes it positive. Thus the signs prescribe the right descents; involutivity gives the same left descents.

For an edge $s\sim t$ incident to a descent, the sufficient inequality $m_s+m_t\leq0$ gives

$$
r_\beta(\alpha_s+\alpha_t)
=\alpha_s+\alpha_t-(m_s+m_t)\beta>0.
$$

This is the geometric terminal test. The construction therefore seeks real roots whose pairings have a specified sign pattern and edge inequalities. In indefinite type, norm two alone is not the real-root witness: the seed is explicitly in a Weyl-group orbit of a simple root, and the iteration is a product of reflections in real roots.

### Why the orbit formula is quadratic

The paper takes $T=r_\gamma r_{\gamma+\delta}$ with $(\gamma,\gamma)=2$ and $(\gamma,\delta)=(\delta,\delta)=0$. If $N=T-1$, composing the reflections gives

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

The null vector belongs to an embedded affine subsystem. It need not be orthogonal to the entire ambient root space; indeed the nonzero pairing $(\beta,\delta)=-a$ produces the quadratic term.

### One further E13 calculation

For $r=3$, the paper’s seed is

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

The positive entries still occur at $I_3=\{0,3,5,7,9,11,12\}$. The exceptional edge sums at $7-8$ and $8-9$ are now $-3$ and $-5$. The general formulas in the paper establish the required inequalities for every nonnegative integer $k$; this single calculation illustrates them.

### Why vanishing needs another argument for covers

A reflection has odd length. If its maximum common descent set has odd size, the only eligible FC non-cover has an even gap, so its coefficient vanishes. Covers have coefficient one and are not eliminated by this reasoning.

The cover lemma considers two possibilities for an FC cover $x$ of $b$. If all descents of $b$ are descents of $x$ on both sides, the maximum-descent lemma forces $x=i(I)$, contradicting the parity of a cover’s length. Otherwise the lifting property forces $x=sb$ or $x=bs$ for a descent $s$ of $b$. The reflection argument then supplies a reduced word containing a braid, contradicting full commutativity.

The assumption that $\beta$ has at least three nonzero coordinates ensures that $s\beta$ is still nonsimple: changing the coordinate at $s$ leaves at least two other nonzero coordinates. This is the place to look for the role of that hypothesis.

The affine family in the paper has a different upper endpoint, an involution that is not a reflection. Its cover exclusion uses the complete FC catalogue and the base-cover computation. Do not transfer the reflection proof to that family without checking its hypotheses.

**Exercise 7.** Derive the formula for $N\rho$ directly from the two reflections, keeping their order. Then derive $N^3=0$ and the orbit formula without assuming that $\delta$ is in the ambient radical.

**Exercise 8.** Recompute $m_2$ and $m_{12}$ from $\beta_{3,1}$ using the branch at node 2. Use the matching printed in the theorem proof to show that $I_3$ is maximum. Explain separately why this gives non-cover vanishing and why the cover lemma is still needed.

## 6 Solution sketches

1. The gaps are $28-4=24$ and $15-4=11$. The first is even, so $\mu=0$ by definition. The second requires $[q^5]P_{x_6,w_6}$, which is one in the displayed $D_6$ polynomial.

2. The inversion counts are $1,4,2,3$. The inverses are respectively $1324,3412,3124,3142$; their descents give the table’s left descents. For $x$, multiplication on the left by $s_1$ gives $e$, outside the star domain, whereas $s_0x=s_0s_1$ is inside. For $w$, $s_1w=s_0s_2s_1$ is inside; $s_0w=s_0s_1s_0s_2s_1=s_1s_0s_1s_2s_1$ has both $s_0$ and $s_1$ as left descents, so is outside.

3. The initial $s$ is a descent of $w$ and hence of $x$. An FC element cannot have noncommuting $s,t$ as simultaneous descents, so $t$ is not a descent of $x$. The condition $d\geq3$ ensures that both possible new gaps, $d$ and $d-2$, stay positive. The gap-one case must already have been handled.

4. The commuting prefix gives the lower bound; support containment and maximum cardinality give the upper bound. The right-cell characterization and its left analogue place $x$ in the two cells of $i(I)$. Their singleton intersection gives equality. Without support containment, a larger commuting factor could use generators outside $K$, beyond the scope of the maximum-cardinality bound.

5. From $e$, left multiplication by $s$ strips back to $e$, and multiplication by $t$ yields $t$. From $t$, multiplication by $s$ yields $st$; multiplication by $t$ yields $e$. From $st$, multiplication by $s$ yields $t$; multiplication by $t$ gives $tst=sts$, which strips on the right to $st$. Right multiplication by an element of $W_J$ preserves the coset $aW_J$. Left multiplication generally changes it.

6. Start from $\{e\}$. At each stage factor any right-terminal $w$ as $av$. Pruning and the induction hypothesis put $v$ in the preceding list; complete coset enumeration supplies $a$, and the exact terminal test accepts $w$. Conversely every accepted product satisfies the definition. Exact equality is needed to recognize repeated representatives and products, and exact descents are needed both for stripping and for the terminal test.

7. Set $u=(\rho,\gamma)$ and $v=(\rho,\delta)$. The first reflection sends $\rho$ to $\rho-(u+v)(\gamma+\delta)$, whose pairing with $\gamma$ is $-u-2v$. Reflecting in $\gamma$ therefore gives $\rho+v\gamma-(u+v)\delta$. The image of $N$ lies in the span of $\gamma,\delta$; $N$ maps that span into the line through $\delta$ and then kills it. Substitution gives the displayed binomial formula.

8. At the branch, $m_2=2(29)-19-26-16=-3$ and $m_{12}=2(16)-29=3$. The matching $(0,1),(2,12),(3,4),(5,6),(7,8),(9,10)$ has six edges, so an independent set uses at most one endpoint of each edge and the one unmatched node: at most seven nodes. The displayed $I_3$ attains seven. Both the reflection length and seven are odd, so the eligible gap is even. A hypothetical FC cover would have coefficient one, which is why the separate cover exclusion remains essential.

## 7 A self-check before moving on

You are ready to read the finite proof without this companion when you can explain why star reduction preserves the coefficient, why maximum common descents force a unique FC endpoint, why parabolic pruning loses no terminal, and exactly which remaining coefficient comes from your thesis.

For the infinite families, add three explanations: why the seed is a real root, why the orbit preserves the terminal inequalities for every parameter, and why no FC cover survives. If an explanation reduces to “the computer checked it,” identify whether the paper actually uses a finite certificate at that point or supplies a symbolic argument. Appendix A records that distinction.
