"""One-time editorial reconstruction from the preserved pre-polish source.

The existing manuscript is edited in place. Mathematical certificates remain
unchanged. This script is retained as an editorial record, not a proof engine.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT/'results/exceptional-leading.tex'
BACKUP = Path(__file__).with_name('source-before-polish.txt')
if not BACKUP.exists():
    BACKUP.write_text(TARGET.read_text())
old = BACKUP.read_text()

def between(a,b):
    return old.split(a,1)[1].split(b,1)[0]

preamble = old.split('\\title{',1)[0]
preamble = preamble.replace('Finite exceptional extensions and infinite generalized E families',
                          'Leading Kazhdan-Lusztig coefficients in type E')
preamble += r'''\newcommand{\supp}{\operatorname{supp}}
\title{Leading Kazhdan--Lusztig coefficients with\newline
fully commutative lower endpoints in type $E$}
\author{}
\date{}
\begin{document}
\maketitle
\begin{abstract}
We study the leading coefficients $\mu(x,w)$ of ordinary equal-parameter
Kazhdan--Lusztig polynomials when $x$ is fully commutative. For finite
Weyl groups of types $E_6$, $E_7$, and $E_8$, we prove
$\mu(x,w)\in\{0,1\}$ for every upper endpoint $w$. The proof combines
star reduction with exhaustive classifications of terminal elements;
parabolic pruning avoids full-group enumeration in type $E_8$.
We also derive a vanishing criterion from maximum commuting descent sets
and construct full-support terminal reflections in every generalized
type $E_{4r+1}$, $r\geq3$. Their leading coefficients with fully
commutative lower endpoint vanish away from covers and, at each fixed
rank, vanish for all sufficiently large values of the family parameter.
Further constructions give complete vanishing for an infinite family
in affine type $\widetilde E_8=E_9$ and non-cover vanishing in type $E_{10}$.
\end{abstract}

\section{Introduction}
Let $(W,S)$ be a Coxeter system, with length $\ell$ and Bruhat order $\leq$.
Write $P_{x,w}(q)$ for its ordinary equal-parameter Kazhdan--Lusztig
polynomials, and define
\[
 \mu(x,w)=\begin{cases}
 [q^{(\ell(w)-\ell(x)-1)/2}]P_{x,w}(q),
       &x<w\text{ and }\ell(w)-\ell(x)\text{ is odd},\\
 0,&\text{otherwise}.
 \end{cases}
\]
A fully commutative element has all its reduced expressions related by
commuting adjacent commuting generators. Denote this set by $W_{\mathrm{FC}}$.
The first result concerns all upper endpoints in the finite exceptional groups.

\begin{theorem}[Finite exceptional types]\label{thm:main}
If $W$ has type $E_6$, $E_7$, or $E_8$, then
\[
 x\in W_{\mathrm{FC}},\quad w\in W
       \quad\Longrightarrow\quad\mu(x,w)\in\{0,1\}.
\]
\end{theorem}
The proof extends Gern's type-$D$ star argument
\cite[Theorem 4.5.11]{Gern}. It reduces to a finite classification of
upper endpoints with no noncommuting initial or final pair in any reduced
expression. We call these elements \emph{terminal}. There are respectively
$1$, $4$, and $6$ noncommuting terminals in $E_6$, $E_7$, and $E_8$.
Each has a unique eligible fully commutative lower endpoint. All resulting
coefficients vanish by parity except for a known pair in a $D_6$ parabolic.
The $E_8$ classification uses exhaustive parabolic pruning rather than
enumeration of its $696,729,600$ elements.

Our second result applies in arbitrarily large ranks. Throughout, the
generalized diagram $E_n$ is the chain $0-1-\cdots-(n-2)$ with node $n-1$
attached to node $2$. Thus $E_9=\widetilde E_8$, and $E_n$ is indefinite
for $n\geq10$. The \emph{support} of an element is the set of generators
occurring in its reduced expressions; full support means $\supp(w)=S$.

\begin{theorem}[A rank-uniform family]\label{thm:uniform}
For every integer $r\geq3$ there is an explicit family of distinct
full-support terminal reflections $b_{r,k}\in W(E_{4r+1})$, $k\geq0$.
For every $x\in W(E_{4r+1})_{\mathrm{FC}}$,
\[
 \mu(x,b_{r,k})\in\{0,1\},\qquad
 x<b_{r,k},\quad \ell(b_{r,k})-\ell(x)>1
       \quad\Longrightarrow\quad\mu(x,b_{r,k})=0.
\]
For each fixed $r$, all these leading coefficients vanish for every
sufficiently large $k$.
\end{theorem}
The construction in Section~\ref{sec:uniform} starts from a staircase
root and iterates an affine-parabolic root action. An even sign change in
$D_{4r}$ proves real-root membership uniformly in $r$. The common descent
set is maximum independent and has odd cardinality. Green's monomial
Temperley--Lieb cell theorem then restricts a potentially nonzero
non-cover coefficient to a single fully commutative lower endpoint,
and parity forces that coefficient to vanish. This proof requires no
Kazhdan--Lusztig polynomial computation.

Green proved the fully commutative lower-endpoint bound in type $A$ and
affine type $A$ \cite{GreenLeading}; Chmutov proved it in finite type $B$
\cite[\S2.3, Theorem 2.3.8]{Chmutov}, hence also for the identical
Coxeter system of type $C$. Green's ADE result with both endpoints fully
commutative \cite[Remark 7.12]{GreenJones} and Jones's result with Deodhar
upper endpoint \cite[Theorem 1.2]{Jones} have different hypotheses.
The structural input for our infinite families is Green's cell theory
in generalized type $E$ \cite{GreenEn}. No claim of historical priority
for the explicit constructions or their consequences is made here.

The infinite-type results concern specified upper families, rather than
arbitrary upper endpoints. They control $\mu$, not all coefficients of
$P_{x,w}$. Section~\ref{sec:applications} gives the complementary affine
$E_8$ and $E_{10}$ families. Section~\ref{sec:transfer} records a
conditional consequence of combinatorial invariance. Exact computational
certificates and reproduction instructions appear in
Appendix~\ref{app:certificates}.

\section{Structural reductions}\label{sec:reductions}
Write
\[
 L(w)=\{s\in S:\ell(sw)<\ell(w)\},\qquad
 R(w)=\{s\in S:\ell(ws)<\ell(w)\}.
\]
An independent vertex set in the Coxeter graph is a set of mutually
commuting simple generators. For such a set $I$, write
$i(I)=\prod_{s\in I}s$; its order is immaterial.
The terminal condition is Gern's ``weakly bad'' condition; a terminal
which is not a commuting product is called ``bad'' in \cite{Gern}.

\subsection{Star reduction}
We use the descent criterion and length-three star transport
\cite[\S2.3 and Theorem 4.2]{KL}, also stated in
\cite[Proposition 1.2.4]{Gern}. Full commutativity is preserved by these
stars in simply laced systems \cite[Proposition 2.10]{Shi}; see also
\cite[Proposition 1.1.29]{Gern}.
'''

reduction = '\\begin{lemma}[Terminal reduction]' + between(
    '\\begin{lemma}[Terminal reduction]',
    '\nThis proof does not import')

structural = r'''
\subsection{Maximum commuting descents}
The next criterion supplies the infinite-family proofs. It does not
require terminality.
'''
structural += '\\begin{lemma}[Maximum common commuting descents]' + between(
    '\\begin{lemma}[Maximum common commuting descents]',
    'For $E_n$ the maximum independent-set size')
structural = structural.replace('let $J$ contain its support', 'let $K$ contain its support')
structural = structural.replace('on $J$', 'on $K$').replace('maximality hypotheses', 'maximum-cardinality hypothesis')
structural += r'''
The independence number of $E_n$ is $\lceil n/2\rceil$. Excluding node $2$
leaves three paths with maximum total $2+\lceil(n-4)/2\rceil$;
including it gives $2+\lceil(n-5)/2\rceil$. Taking the larger proves the claim.

\subsection{An ambient root action}\label{subsec:rootaction}
Normalize simple roots by $(\alpha_j,\alpha_j)=2$ and
$(\alpha_s,\alpha_t)=-1$ for adjacent nodes. For a real root $\beta$,
its reflection acts by $r_\beta(\rho)=\rho-(\rho,\beta)\beta$.
Root inequalities refer to positivity in the simple-root basis.

\begin{lemma}\label{lem:rootaction}
Suppose that $\gamma$ and $\gamma+\delta$ are real roots, with
$(\gamma,\gamma)=2$ and $(\gamma,\delta)=(\delta,\delta)=0$.
For $T=r_\gamma r_{\gamma+\delta}$,
\[
 (T-1)\rho=(\rho,\delta)\gamma-
       ((\rho,\gamma)+(\rho,\delta))\delta,\qquad (T-1)^3=0.
\]
If $(\beta,\delta)=-a$ and $(\beta,\gamma)=-c$, then
\begin{equation}\label{eq:rootiterate}
 T^k\beta=\beta-ak\gamma+(ak^2+ck)\delta\qquad(k\geq0).
\end{equation}
\end{lemma}
\begin{proof}
The first identity follows by multiplying the two reflection formulas.
Set $N=T-1$. Then $N\delta=0$ and $N\gamma=-2\delta$, so $N^3=0$.
Moreover $N\beta=-a\gamma+(a+c)\delta$ and $N^2\beta=2a\delta$.
The binomial expansion of $(1+N)^k$ gives \eqref{eq:rootiterate}.
\end{proof}
The lemma applies on the whole ambient root space. In the indefinite
constructions below, $\delta$ is null in an affine parabolic but is not
in the radical of the ambient form.
'''

u = between('\\section{A rank-uniform full-support construction}',
            '\\section{Use of the new invariance input}')
notation = u.split('\\begin{theorem}',1)[0]
notation = notation.replace('\\label{sec:uniform}','')
notation = notation.replace('The $E_{13}$ seed extends to every rank $n=4r+1$ with $r\\geq3$.\nThis statement uses the whole diagram in each rank; it is not obtained\nby embedding a fixed smaller upper endpoint. ', '')
notation = notation.replace('\\tag{***}', '\\label{eq:uniformroots}')
notation = notation.replace('\\[\n \\beta_{r,k}', '\\begin{equation}\n \\beta_{r,k}')
notation = notation.replace('\\label{eq:uniformroots}\n\\]', '\\label{eq:uniformroots}\n\\end{equation}')
uniform = '\\section{The rank-uniform construction}\\label{sec:uniform}\n' + notation
up = u.split('\\begin{proof}',1)[1].split('\\end{proof}',1)[0]
seed = up.split('The roots $\\gamma$',1)[0]
seed = seed.replace('We first prove real-root membership uniformly, rather than inferring it\nfrom the norm. ', '')
uniform += r'''
\subsection{The staircase seed}
\begin{lemma}\label{lem:uniformseed}
For every $r\geq3$, the vector $\beta_r$ is a real root of $E_{4r+1}$.
\end{lemma}
\begin{proof}
''' + seed + '\\end{proof}\n'
middle = 'The roots $\\gamma$' + up.split('The roots $\\gamma$',1)[1].split(
    'The independent sets $I_{10}$',1)[0]
middle = middle.split('The independent set',1)[0] if 'The independent set' in middle else middle
# The original final proof starts at the independence-number assertion.
middle = middle.split('The graph has independence number',1)[0]
middle = middle.replace('the earlier operator identity\nholds on the whole ambient space',
                        'Lemma~\\ref{lem:rootaction} applies')
middle = middle.replace('The nilpotent identity $(T-1)^3=0$ consequently gives',
                        'Equation~\\eqref{eq:rootiterate} gives')
middle = middle.replace('in (***)', 'in \\eqref{eq:uniformroots}')
middle = middle.replace('The reflection formula now gives both descent sets $I_r$ and\nterminality exactly as in the fixed-rank proof.', r'''For $s\in I_r$ and $t\sim s$,
\[
 b_{r,k}(\alpha_s+\alpha_t)
 =\alpha_s+\alpha_t-(m_s+m_t)\beta_{r,k}>0.
\]
Thus deleting a right descent creates no adjacent right descent.
The reflection is an involution, so the same holds on the left.''')
uniform += r'''
\subsection{Translated roots and descent sets}
\begin{lemma}\label{lem:uniformdescents}
For all $r\geq3$ and $k\geq0$, $\beta_{r,k}=T^k\beta_r$ is a real root
with every coordinate strictly positive. Its reflection is terminal and
has both descent sets equal to $I_r$.
\end{lemma}
\begin{proof}
''' + middle + '\\end{proof}\n'
tail = up.split('The graph has independence number',1)[1]
tail = 'The graph has independence number' + tail
tail = tail.replace('The graph has independence number $2r+1$, attained by $I_r$.',r'''The set $I_r$ is independent of size $2r+1$.
The matching
\[
 (0,1),(2,4r),(3,4),(5,6),\ldots,(4r-3,4r-2)
\]
has $2r$ edges and leaves one vertex unmatched, proving that $I_r$ is maximum.''')
tail = tail.replace('without computing any KL polynomial', '')
uniform += '\n\\begin{proof}[Proof of Theorem~\\ref{thm:uniform}]\n' + tail + '\\end{proof}\n'
uniform += r'''
At $r=3$ the seed is
\[
 \beta_3=(2,3,5,5,4,4,3,3,2,2,1,1,3),
\]
giving the $E_{13}$ family. The restriction $r\geq3$ is structural:
at $r=2$, $(\alpha_0,\beta_2)=0$ and the maximum odd descent set is lost.
No effective threshold for the eventual vanishing in
Theorem~\ref{thm:uniform} is asserted.

\section{Finite exceptional types}\label{sec:finite}
\subsection{The terminal classification}
A fully commutative $x$ is \emph{eligible} for a terminal $b$ if
\[
 x\leq b,\qquad L(b)\subseteq L(x),\qquad R(b)\subseteq R(x).
\]
For non-covers, every ineligible pair has $\mu(x,b)=0$ by the descent
criterion. Words in the finite types are written as strings of node labels.

\begin{proposition}[Computational classification]\label{prop:classification}
The noncommuting terminals in $E_6$, $E_7$, and $E_8$ are exactly those
in Table~\ref{tab:terminals}. Each has the unique eligible lower endpoint
shown, and both its descent sets are the support of that endpoint.
The remaining terminals are commuting products: $22$, $36$, and $58$,
respectively, including the identity.
\end{proposition}
\begin{table}[ht]
\centering\small
\begin{tabular}{@{}clrrr@{}}
\toprule
Type & Terminal $b$ & $\ell(b)$ & Eligible $x$ & $\ell(b)-\ell(x)$\\
\midrule
$E_6$ & $1325213$ & 7 & $135$ & 4\\
$E_7$ & $1326213$ & 7 & $136$ & 4\\
$E_7$ & $13256213$ & 8 & $1356$ & 4\\
$E_7$ & $132543621324356$ & 15 & $1356$ & 11\\
$E_7$ & $1325436210321432543621324356$ & 28 & $1356$ & 24\\
$E_8$ & $1327213$ & 7 & $137$ & 4\\
$E_8$ & $13257213$ & 8 & $1357$ & 4\\
$E_8$ & $61327213$ & 8 & $1367$ & 4\\
$E_8$ & $132543721324357$ & 15 & $1357$ & 11\\
$E_8$ & $1325437210321432543721324357$ & 28 & $1357$ & 24\\
$E_8$ & $b_{8,6}$ & 50 & $1357$ & 46\\
\bottomrule
\end{tabular}
\caption{All noncommuting terminals and their eligible lower endpoints.}
\label{tab:terminals}
\end{table}
Here $b_{8,6}$ is the concatenation, without cancellation, of
\[
 (7534231270123456210321432)(5437210321432543721324357).
\]

\subsection{Exhaustiveness of the classification}
All computations use the faithful integral geometric representation.
For columns $a_j=w(\alpha_j)$, right multiplication by $s_i$ replaces
\[
 a_i\longmapsto-a_i,\qquad a_j\longmapsto a_j+a_i\quad(j\sim i),
\]
and fixes the other columns. Signs determine descents. The terminal test is
\begin{equation}\label{eq:terminaltest}
 s\in R(w)\quad\Longrightarrow\quad
 w(\alpha_s+\alpha_t)>0\quad\text{for every }t\sim s,
\end{equation}
together with its counterpart for $w^{-1}$. Indeed, a noncommuting suffix
$ts$ is equivalent to $s\in R(w)$ and $t\in R(ws)$.
In $E_6$ and $E_7$, breadth-first multiplication by all generators until
closure enumerates the complete group and gives minimal word lengths.

For $E_8$, write the unique length-additive parabolic decomposition
$w=av$, with $a\in W^J$ and $v\in W_J$, where $W^J$ consists of the
minimal right coset representatives. Call an element \emph{right-terminal}
if it has no reduced expression ending in a noncommuting pair.

\begin{lemma}[Parabolic pruning]\label{lem:pruning}
If $w=av$ is right-terminal, then $v$ is right-terminal in $W_J$.
\end{lemma}
\begin{proof}
A reduced expression for $v$ ending in adjacent generators would,
after prefixing a reduced expression for $a$, give such an expression
for $w$. Length additivity rules this out.
\end{proof}
Consequently one may enumerate all products $av$ with $v$ right-terminal,
apply \eqref{eq:terminaltest} in the ambient group, and then test inverses
to select the two-sided terminals. The recursive implementation follows
\[
 E_8\supset E_7\supset E_6\supset D_5\supset D_4\supset A_3
 \supset A_2\supset A_1\supset\{e\}.
\]
At each step all minimal coset representatives are obtained by closure
under left multiplication, reducing to right-$J$-minimal form by stripping
right $J$-descents. This reaches every coset because the left action on
right cosets is transitive. Starting with the identity, induction using
Lemma~\ref{lem:pruning} proves completeness at every stage.
The final $E_7$ parabolic has $576$ right-terminals and $240$ cosets.
Testing their products gives $2,160$ right-terminals in $E_8$; the inverse
test leaves exactly the $64$ terminals in
Proposition~\ref{prop:classification}.

Full commutativity is tested by the recurrence
\begin{equation}\label{eq:FCrecurrence}
 w\in W_{\mathrm{FC}}\quad\Longleftrightarrow\quad
 R(w)\text{ is commuting and }ws\in W_{\mathrm{FC}}
      \text{ for every }s\in R(w).
\end{equation}
Matsumoto's theorem justifies this test: a braid before the last letter
survives in a shorter prefix, while a length-three braid ending at the
last letter gives adjacent right descents. For $E_8$ it suffices to close
the FC set under all ascents passing \eqref{eq:FCrecurrence}. Every FC
element has FC weak predecessors, so a minimal omitted element would
contradict closure. This gives the complete $10,846$-element catalogue.
For every terminal, a Bruhat comparison using the lifting property and
both descent inclusions exhausts all eligible endpoints in that catalogue.
The same procedure is applied in $E_6$ and $E_7$. Independent implementations,
counts, and certificate locations are given in
Appendix~\ref{app:certificates}.

\begin{proof}[Proof of Theorem~\ref{thm:main}]
By Lemma~\ref{lem:reduction}, it suffices to treat terminal upper endpoints.
For a commuting product $b$, every $x<b$ is a proper subproduct. The
descent criterion gives zero on non-covers and one on covers.
For a noncommuting terminal, Table~\ref{tab:terminals} and parity leave
only the length-$15$ rows. These are the same pair in a $D_6$ parabolic.
Mapping Gern's labels $1,2,3,4,5,6$ to $1,n-1,2,3,4,5$ identifies the
eligible pair with
\[
 x_6=(-1,-2,4,3,6,5),\qquad w_6=(-1,-6,3,-4,5,-2).
\]
Bruhat order and ordinary KL polynomials agree with those in the standard
parabolic. Gern's Lemma~4.5.5 \cite{Gern} gives $\mu(x_6,w_6)=1$.
All other non-cover coefficients are zero, proving the theorem.
\end{proof}

\section{Complementary fixed-rank families}\label{sec:applications}
'''

affine = between('\\section{An unbounded full-support affine-$E_8$ family}',
                 '\\section{Full-support infinite families')
affine = affine.split('The all-parameter proof and an independent verifier',1)[0]
affine = '\\subsection{Complete vanishing in affine type $\\widetilde E_8$}\n' + affine
matrix = '\\[\n B_0=' + affine.split('\\[\n B_0=',1)[1].split('\\]\n',1)[0] + '\\]\n'
affine = affine.replace(matrix, 'The columns are recorded in Appendix~\\ref{app:affinedata}.\n')
affine = affine.replace('The explicit word gives\n$b_0^2=1$ and the following matrix, whose columns are $b_0(\\alpha_j)$:',
                         'Multiplying the displayed word gives $b_0^2=1$.\nWrite $B_0$ for the matrix with columns $b_0(\\alpha_j)$.')
affine = affine.replace(',&x_I&=i(I)', '')
affine = affine.replace('In particular these are infinitely many distinct upper endpoints in one\nfixed infinite Coxeter group. They are not a whole-type classification.\n', '')
affine = affine.replace('\\tag{*}', '\\label{eq:affinecolumns}')
affine = affine.replace('\\[\n b_k(\\alpha_j)', '\\begin{equation}\n b_k(\\alpha_j)')
affine = affine.replace('\\label{eq:affinecolumns}\n\\]', '\\label{eq:affinecolumns}\n\\end{equation}')
affine = affine.replace('Consequently (*)', 'Consequently \\eqref{eq:affinecolumns}')
affine = affine.replace('These finite arithmetic checks\nprove', 'These checks prove')
affine = affine.replace('All $240$ rows are saved and independently reproduced;\nthis is an all-parameter inversion count, not a fitted length formula.',
                         'The complete root-string data are included in the certificate. ')
affine = affine.replace('The independent complete FC\ncatalogue has maximum length $44$',
                         'The complete FC catalogue has maximum length $44$')
affine = affine.replace('This last finite\ncheck eliminates the cover exception and proves complete vanishing.',
                         'This eliminates the cover exception. If $b_k$ were FC,\nLemma~\\ref{lem:maxI} would force $b_k=i(I)$, contradicting\n$\\ell(b_k)=27+92k>5$.')
affine = affine.replace('These finite arithmetic checks', 'These direct checks')

e10 = r'''
\subsection{An indefinite $E_{10}$ family}
Embed the affine subdiagram on $\{0,\ldots,7,9\}$ and take
\[
 \delta=(2,4,6,5,4,3,2,1,0,3),\qquad
 \gamma=(1,2,3,2,2,1,1,0,0,1).
\]
Here $\gamma$ differs from the root used in the preceding affine family.
Set
\[
 \beta_0=(3,7,10,9,7,6,4,3,1,6),\qquad I=\{1,3,5,7,9\},
\]
and define
\begin{equation}\label{eq:e10roots}
 \beta_k=\beta_0-k\gamma+(k^2+3k)\delta,\qquad b_k=r_{\beta_k}.
\end{equation}
\begin{proposition}\label{thm:indefinite}
The reflections in \eqref{eq:e10roots} are distinct full-support terminal
elements of $W(E_{10})$, with $L(b_k)=R(b_k)=I$. For every fully
commutative $x$, $\mu(x,b_k)\in\{0,1\}$, and every non-cover coefficient
vanishes.
\end{proposition}
\begin{proof}
The lowering word in Appendix~\ref{app:affinedata} reduces $\beta_0$
to $\alpha_9$, proving real-root membership. The roots $\gamma$ and
$\gamma+\delta$ are real roots of the embedded affine parabolic, with
$(\beta_0,\gamma)=-3$ and $(\beta_0,\delta)=-1$.
Lemma~\ref{lem:rootaction} gives $\beta_k=T^k\beta_0$.
The binomial increments $4\delta-\gamma$ and $2\delta$ are nonnegative,
so every coordinate of $\beta_k$ is strictly positive.
Its simple-root pairings are
\begin{align*}
 ((\alpha_j,\beta_k))_{j=0}^9
 ={}&(-1,1,-2-k,1+k,-1-k,1+k,-1-k,1+k,\\
    &\hspace{30mm}-1-k^2-3k,2+k).
\end{align*}
The positive entries occur exactly on $I$, and their sums with adjacent
entries are nonpositive. The reflection formula and
\eqref{eq:terminaltest} give both descent sets and terminality.
The constant pairing $(\alpha_1,\beta_k)=1$ makes column $1$ change
every coordinate row, proving full support. The coordinate
$\beta_{k,0}=2k^2+5k+3$ proves distinctness.
The set $I$ is maximum independent of size $5$. Reflection parity and
Proposition~\ref{prop:parity} give the coefficient assertions.
As in Theorem~\ref{thm:uniform}, the reflections are not FC because their
$(-1)$-eigenspaces have dimension one, whereas that of $i(I)$ has dimension five.
\end{proof}
The base reflection has length $101$. No formula for the lengths of the
indefinite families is needed. Standard-parabolic embeddings transport
the $E_9$ and $E_{10}$ families to larger generalized $E_n$ groups;
their supports remain those of the smaller parabolics.

\section{Interval transfer and scope}\label{sec:transfer}
Assume the combinatorial-invariance theorem of \cite[Theorem 1.1]{CI}:
isomorphic closed Bruhat intervals have the same ordinary equal-parameter
KL polynomial. It gives the following consequence of Theorem~\ref{thm:main}.

\begin{corollary}[Conditional interval transfer]
If a closed Bruhat interval $[u,v]$ in a Coxeter group is isomorphic as
a poset to $[x,w]$ in finite $E_6$, $E_7$, or $E_8$, with $x$ fully
commutative, then $\mu(u,v)\in\{0,1\}$.
\end{corollary}
\begin{proof}
Combinatorial invariance gives $P_{u,v}=P_{x,w}$. The interval isomorphism
preserves relative rank, hence the exponent defining $\mu$.
Apply Theorem~\ref{thm:main}.
\end{proof}
The same argument transfers the bounds of Theorem~\ref{thm:uniform} and
the fixed-rank families to intervals admitting those models. In particular,
an affine-family model forces $\mu(u,v)=0$. The target lower endpoint
need not be fully commutative. This is a conditional transfer statement;
no further family of cross-type interval isomorphisms is asserted.

Vanishing of $\mu$ does not imply $P=1$. In $D_4$, let $x$ be the product
of the three leaves and let $c$ be the center. The element $b=xcx$ has
maximum common commuting descents and an even gap $\ell(b)-\ell(x)=4$,
but $P_{x,b}=1+2q$ \cite[\S3, preprint Corollary 3]{Mongelli}.
Thus the distinction between the leading coefficient and the full
polynomial is essential, even for the descent mechanism used here.

\appendix
\section{Computational certificates}\label{app:certificates}
The finite classification and the affine cover check are computational
parts of the proofs. All arithmetic is exact. For $E_6$ and $E_7$, one
implementation stores simple-root images as root indices, and an
independent implementation stores integer columns; neither imports the other.
For $E_8$, two different parabolic decompositions certify the same terminal
list. The $E_7$ method uses $138,240$ candidates at its final stage and
$143,660$ across the entire recursive chain. The independent $D_7$ method
uses $280$ right-terminals and $2,160$ cosets, giving $604,800$ candidates.

The latter cosets are obtained from the orbit of the integral fundamental
weight $(4,7,10,8,6,4,2,5)$. Its simple-root pairings are $1$ at node $0$
and zero elsewhere. Distinct orbit images give distinct cosets; the
orbit has the full index $2,160$, and the stored representatives have
no right descent in $D_7$. These checks establish completeness and minimality.

\begin{center}
\begin{tabular}{@{}lrrr@{}}
\toprule
Check & $E_6$ & $E_7$ & $E_8$\\
\midrule
Roots & 72 & 126 & 240\\
FC elements & 662 & 2,670 & 10,846\\
Noncommuting terminals & 1 & 4 & 6\\
Commuting terminals & 22 & 36 & 58\\
FC star moves checked & 3,620 & 16,776 & 77,476\\
\bottomrule
\end{tabular}
\end{center}
The full $E_6$ and $E_7$ searches contain $51,840$ and $2,903,040$
elements, with maximum lengths $36$ and $63$. Their length distributions
agree with the Weyl-group Poincar\'e products. The FC counts agree with
\cite[\S5.1]{BJN}. These are consistency checks; exhaustiveness follows
from the closure and pruning arguments in Section~\ref{sec:finite}.
An independent $E_8$ verifier checks the terminal words across both
decompositions, lengths by all $120$ positive-root inversions, both terminal
tests, complete FC closure, and every eligible lower endpoint. It also
checks the $D_7$ parabolic/coset Poincar\'e product against the $E_8$
exponents $1,7,11,13,17,19,23,29$.

The affine $E_9$ FC catalogue has $44,199$ elements and maximum length
$44$. Its closure is certified by \eqref{eq:FCrecurrence}. Deleting one
letter from the length-$27$ base word gives exactly $21$ distinct Bruhat
covers of length $26$, none fully commutative. Independent root matrices
also verify all $240$ affine inversion-string records used in
Theorem~\ref{thm:affine}.

The certificate archive is organized as follows:
\begin{center}\small
\begin{tabular}{@{}ll@{}}
\toprule
Result & Programs and certificates\\
\midrule
$E_6,E_7$ & \path{research/verify_exceptional.py}\\
 & \path{results/e6-independent-certificate.json}\\
 & \path{results/e7-independent-certificate.json}\\
$E_8$ & \path{research/en_e8/reproduce.py}\\
 & \path{research/en_independent/verify_e8.py}\\
 & \path{results/e8-independent-audit.json}\\
Affine $E_9$ & \path{research/en_affine_referee/verify_families.py}\\
 & \path{research/en_affine_referee/fc-cover-certificate.json}\\
Uniform family & \path{research/en_uniform/construction.py}\\
 & \path{research/en_uniform/referee_verify.py}\\
 & \path{research/en_uniform/referee-certificate.json}\\
$E_{10}$ & \path{research/en_affine_referee/verify_indefinite_e10.py}\\
 & \path{research/en_affine_referee/e10-certificate.json}\\
\bottomrule
\end{tabular}
\end{center}
The uniform root proof is symbolic. Arithmetic checks in finitely many
ranks are supplementary, and do not replace the all-rank argument.
The separately verified $E_{13}$ certificate is retained in
\path{research/en_affine_referee/e13-certificate.json}.
The commands \path{python3 research/reproduce_en.py} and
\path{python3 research/reproduce_uniform.py} reproduce the proof
certificates using Python's standard library and, for the finite searches,
an installed C++17 compiler. Source hashes, logs, and complete code are
included in the archive. The focused literature audit is
\path{research/en_uniform/literature.txt}.

\section{Seed data}\label{app:affinedata}
The matrix of the affine base word, with columns $b_0(\alpha_j)$, is
'''

witness = r'''
For the $E_{10}$ seed, successive simple reflections in the following
order reduce $\beta_0$ to $\alpha_9$:
\begin{align*}
 &1,3,5,4,7,6,5,9,2,1,0,3,2,1,4,3,2,\\
 &5,6,7,8,9,2,1,0,3,2,1,4,3,2,5,4,3,6,\\
 &5,4,7,6,5,9,2,1,0,3,2,1,4,3,2.
\end{align*}
Each intermediate vector is nonnegative and its height decreases.
For the uniform finite-$E_8$ root $\gamma$, the fixed sequence
\[
 2,1,0,4,3,2,1,6,5,4,3,2
\]
reduces it to the branch simple root. For $\gamma+\delta$, append
\begin{align*}
 &b,2,1,0,3,2,1,4,3,2,5,4,3,6,5,4,7,6,5,\\
 &b,2,1,0,3,2,1,4,3,2,
\end{align*}
where $b$ denotes the branch label ($4r$ or $9$). These fixed witnesses
also certify the affine real-root input to Lemma~\ref{lem:rootaction}.
'''

bib = '\\begin{thebibliography}{99}' + old.split('\\begin{thebibliography}{99}',1)[1]
for key in ('duCloux','Woo','BG'):
    bib = re.sub(r'\\bibitem\{'+key+r'\}.*?(?=\\bibitem|\\end\{thebibliography\})',
                 '',bib,flags=re.S)
text = preamble + reduction + structural + uniform + affine + e10 + matrix + witness + bib
TARGET.write_text(text)
print(f'Revised {TARGET}: {len(text.splitlines())} lines; original preserved at {BACKUP}.')
