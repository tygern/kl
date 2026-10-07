% LaTeX fragment for inclusion in the parent research note; not a standalone document.
\section{An exhaustive exceptional-type extension}
\label{sec:exceptional-direct}
We use the numbering in which $0,1,\ldots,r-2$ form a chain
and $r-1$ is attached to $2$.  Thus the edges are
$\{i,i+1\}$ for $0\leq i<r-2$, together with $\{2,r-1\}$.
Words below are strings of these labels, interpreted from left to right.

\begin{theorem}[Computer-assisted finite-type extension]
For $W$ of type $E_6$ or $E_7$, if $x$ is fully commutative and $w\in W$
is arbitrary, then $\mu(x,w)\in\{0,1\}$.
\end{theorem}

The proof uses exhaustive enumeration only to classify a very small set of
terminal elements and their possible lower endpoints.  It does not compute
all Kazhdan--Lusztig polynomials, and does not invoke combinatorial invariance.
Historical priority for the complete exceptional-type assertion is not certified
by the literature search accompanying this computation.

\begin{lemma}[Terminal reduction]
Let $W$ be simply laced.  Suppose $\mu(x,b)\in\{0,1\}$ whenever $x$ is fully
commutative and $b$ has no reduced expression beginning or ending with a pair
of noncommuting generators.  Then the same conclusion holds for arbitrary
upper endpoints in $W$.
\end{lemma}
\begin{proof}
Induct on $\ell(w)$.  There is nothing to prove if $x\not<w$, if the length
difference is even, or if the length difference is one: the corresponding
values are respectively $0$, $0$, and $1$.  We may therefore assume
$d=\ell(w)-\ell(x)\geq3$ is odd.
If $L(w)\not\subseteq L(x)$ or $R(w)\not\subseteq R(x)$, the usual descent
criterion gives $\mu(x,w)=0$; the exceptional possibility in that criterion
is a cover, which has already been excluded.  Hence assume both inclusions.
Since the descent sets of a fully commutative element are pairwise commuting,
so are $L(w)$ and $R(w)$.

If $w$ is terminal we use the hypothesis.  Otherwise some reduced expression
begins or ends with noncommuting generators.  Consider the first case,
$w=stv$ reduced, where $m(s,t)=3$.  We have $s\in L(w)$, $t\notin L(w)$,
and the left $\{s,t\}$-star operation satisfies ${}^*w=sw$, of length
$\ell(w)-1$.  The descent inclusions give $s\in L(x)$; full commutativity
then gives $t\notin L(x)$.  Consequently ${}^*x$ is defined and is fully
commutative.  The Kazhdan--Lusztig star-transport identity gives
\[
 \mu(x,w)=\widetilde\mu({}^*x,{}^*w),
\]
where $\widetilde\mu$ is the symmetric version, equal to zero for incomparable
pairs.  The transformed length difference is $d$ or $d-2$, hence remains
positive.  Thus either the new endpoints are incomparable and the value is
zero, or ${}^*x<{}^*w$ and the induction hypothesis applies to the strictly
shorter upper endpoint ${}^*w$.  The right-handed case is identical.
\end{proof}

\paragraph{Exact terminal classification.}
Write $[i_1\cdots i_k]=s_{i_1}\cdots s_{i_k}$.  In type $E_6$ there are $22$
terminal products of commuting generators and exactly one other terminal:
\[
 b=[1325213],\qquad \ell(b)=7,\qquad L(b)=R(b)=\{1,3,5\}.
\]
The only fully commutative $x\leq b$ satisfying
$L(b)\subseteq L(x)$ and $R(b)\subseteq R(x)$ is $x=[135]$, of length $3$.
The relevant difference is $4$, so its $\mu$-value is zero.

In type $E_7$ there are $36$ terminal products of commuting generators and
exactly the following four additional terminals:
\[
\begin{array}{c|l|c|c}
 & \text{reduced word for } b & \ell(b)&\text{unique eligible }x\\\hline
b_1&1326213&7&[136]\\
b_2&13256213&8&[1356]\\
b_3&132543621324356&15&[1356]\\
b_4&1325436210321432543621324356&28&[1356]
\end{array}
\]
Here ``eligible'' means fully commutative, below $b$, and containing both
descent sets of $b$ in the corresponding descent sets of $x$.
For $b_1$, both descent sets are $\{1,3,6\}$; for $b_2,b_3,b_4$, both
are $\{1,3,5,6\}$.  The four length differences are respectively
$4,4,11,24$.  Parity disposes of the first, second, and fourth cases.
The third is in the standard $D_6$ parabolic on $\{1,2,3,4,5,6\}$, so
Gern's Theorem~4.5.11 applies.  More precisely, mapping Gern's generators
$1,2,3,4,5,6$ to $1,6,2,3,4,5$ identifies this pair with
\[
 x=(-1,-2,4,3,6,5),\qquad b_3=(-1,-6,3,-4,5,-2),
\]
for which Gern's Lemma~4.5.5 gives $\mu(x,b_3)=1$.
Products of commuting generators give only $\mu$-values $0$ and $1$ by
the descent criterion.  Terminal reduction now proves the theorem.

\paragraph{Why the finite computation is exhaustive.}
The attached program constructs the roots in integral simple-root coordinates
by closing the simple roots under the reflections
\[
 (s_i\alpha)_i=-\alpha_i+\sum_{j\sim i}\alpha_j,
 \qquad (s_i\alpha)_j=\alpha_j\quad(j\ne i).
\]
An element is represented faithfully by its images of the simple roots.
Breadth-first multiplication by every simple reflection enumerates the whole
group and determines Coxeter length.  The programs obtain respectively
$72,126$ roots, $51\,840,2\,903\,040$ elements, and longest lengths $36,63$.
A right descent is tested independently by the sign of the corresponding
simple-root image, and is checked against the computed length.

The terminal test is exactly
\[
 \forall s\in R(w),\quad R(ws)\cap\{t:t\sim s\}=\varnothing,
\]
and its left-handed analogue, obtained using inversion.  This is equivalent
to the stated reduced-expression definition, because a reduced suffix $ts$
is equivalent to $s\in R(w)$ and $t\in R(ws)$.
Full commutativity is computed inductively: $w$ is fully commutative if and
only if all $ws$, $s\in R(w)$, are fully commutative and $R(w)$ is commuting.
Indeed, a forbidden braid occurring before the last letter is detected in a
shorter prefix; a forbidden braid ending at the last letter gives two adjacent
right descents.  This also proves the converse.
The resulting fully commutative counts are $662$ and $2670$.
Finally, the lifting-property comparison checks each eligible bottom against
each of the finitely many terminals.  Independent integral-matrix enumerations
in the parent project reproduce the same lists.

\paragraph{Scope.}
This is a direct, computer-assisted extension to two entire exceptional Weyl
groups.  It is not a construction of new $D$-model Bruhat intervals.  In
particular, star operations preserve the relevant $\mu$-coefficient; no claim
that they preserve complete interval posets or full Kazhdan--Lusztig
polynomials is made.  The exceptional terminal $b_4$ has full $E_7$ support,
so the terminal classification goes beyond listing classical parabolic tops,
even though parity handles its remaining $\mu$-value.
