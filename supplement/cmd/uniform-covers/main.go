// Command uniform-covers performs finite checks of the ingredients of the
// manuscript's Lemma on covers of terminal reflections, for the rank-uniform
// family b_{r,k} = r_{beta_{r,k}} in W(E_{4r+1}), and the observed length
// formula l(b_{r,k}) = 8r^2 + 3 + 116k.
//
// For each (r,k) checked:
//   - beta_{r,k} has norm 2, b = r_beta is a terminal involution of odd length
//     8r^2 + 3 + 116k with L(b) = R(b) = I_r, full support, not fully
//     commutative (Stembridge's heap criterion on a reduced word);
//   - for every s in I_r the element b s has L(b s) = I_r, R(b s) = I_r minus {s},
//     length l(b) - 1, and is not fully commutative; s b is its inverse;
//   - s b s = r_{s beta} with s beta a positive non-simple root, and
//     (alpha_s, s beta) = -(alpha_s, beta) < 0 (the pairing used in Lemma A);
//   - for selected (r,k): every Bruhat cover of b (one letter deleted from a
//     reduced word, keeping length l(b) - 1; exhaustive by the strong exchange
//     condition) is not fully commutative.
//
// Ports research/review_checks/uniform_family_checks.py.  Standard library
// only; it imports no other package of this module (self-contained review
// program).  Elements are integer matrices stored as columns (column j =
// w(alpha_j) in the simple-root basis), exact integer arithmetic with overflow
// guards.  E_n labelling: chain 0-1-...-(n-2), node n-1 attached to node 2.
// A string of labels denotes the product of the simple reflections in the
// order written.
//
// Usage: uniform-covers [-full]
// Default: r = 3..8, k = 0..3 and covers for (r,k) in {(3,0),(3,1),(4,0)}.
// -full adds r = 9, 10, 12 and covers for (3,2),(4,1),(5,0),(6,0),(7,0).
// Writes results/uniform-family-certificate.json (relative to the current
// directory).
package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

func addChk(a, b int64) int64 {
	c := a + b
	if (b > 0 && c < a) || (b < 0 && c > a) {
		panic("int64 overflow in addition")
	}
	return c
}

func mulChk(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	c := a * b
	if c/b != a || (a == -1 && b == math.MinInt64) || (b == -1 && a == math.MinInt64) {
		panic("int64 overflow in multiplication")
	}
	return c
}

// Elem is an n x n integer matrix stored column-major: entry (i,j) = d[j*n+i].
type Elem []int64

// En is the Weyl group of type E_n (simply laced Dynkin diagram).
type En struct {
	n      int
	adj    [][]int
	adjset []map[int]bool
	e      Elem
}

func newEn(n int) *En {
	g := &En{n: n}
	g.adj = make([][]int, n)
	for i := 0; i < n-2; i++ {
		g.adj[i] = append(g.adj[i], i+1)
		g.adj[i+1] = append(g.adj[i+1], i)
	}
	g.adj[n-1] = append(g.adj[n-1], 2)
	g.adj[2] = append(g.adj[2], n-1)
	g.adjset = make([]map[int]bool, n)
	for i := range g.adj {
		g.adjset[i] = map[int]bool{}
		for _, j := range g.adj[i] {
			g.adjset[i][j] = true
		}
	}
	g.e = make(Elem, n*n)
	for j := 0; j < n; j++ {
		g.e[j*n+j] = 1
	}
	return g
}

// pair computes (v, alpha_i) for the symmetric form with (alpha_i, alpha_i) = 2.
func (g *En) pair(v []int64, i int) int64 {
	s := mulChk(2, v[i])
	for _, j := range g.adj[i] {
		s = addChk(s, -v[j])
	}
	return s
}

func (g *En) form(v, w []int64) int64 {
	var s int64
	for i := 0; i < g.n; i++ {
		s = addChk(s, mulChk(v[i], g.pair(w, i)))
	}
	return s
}

func (g *En) reflVec(v []int64, i int) []int64 {
	c := g.pair(v, i)
	out := make([]int64, len(v))
	copy(out, v)
	if c == 0 {
		return out
	}
	out[i] = addChk(out[i], -c)
	return out
}

func (g *En) col(w Elem, j int) []int64 { return w[j*g.n : (j+1)*g.n] }

func (g *En) clone(w Elem) Elem {
	out := make(Elem, len(w))
	copy(out, w)
	return out
}

// rmulInPlace replaces w by w * s_i.
func (g *En) rmulInPlace(w Elem, i int) {
	n := g.n
	ai := g.col(w, i)
	for _, j := range g.adj[i] {
		cj := w[j*n : (j+1)*n]
		for k := 0; k < n; k++ {
			cj[k] = addChk(cj[k], ai[k])
		}
	}
	for k := 0; k < n; k++ {
		ai[k] = -ai[k]
	}
}

func (g *En) rmul(w Elem, i int) Elem {
	out := g.clone(w)
	g.rmulInPlace(out, i)
	return out
}

func (g *En) lmul(i int, w Elem) Elem {
	out := make(Elem, len(w))
	for j := 0; j < g.n; j++ {
		copy(out[j*g.n:(j+1)*g.n], g.reflVec(g.col(w, j), i))
	}
	return out
}

func (g *En) mul(u, v Elem) Elem {
	n := g.n
	out := make(Elem, n*n)
	for j := 0; j < n; j++ {
		col := out[j*n : (j+1)*n]
		for i := 0; i < n; i++ {
			c := v[j*n+i]
			if c != 0 {
				ui := u[i*n : (i+1)*n]
				for k := 0; k < n; k++ {
					col[k] = addChk(col[k], mulChk(c, ui[k]))
				}
			}
		}
	}
	return out
}

func negCol(col []int64) bool {
	any := false
	for _, x := range col {
		if x < 0 {
			any = true
		}
		if x > 0 {
			return false
		}
	}
	return any
}

func posCol(col []int64) bool {
	any := false
	for _, x := range col {
		if x > 0 {
			any = true
		}
		if x < 0 {
			return false
		}
	}
	return any
}

// R is the right descent set (increasing list).
func (g *En) R(w Elem) []int {
	out := []int{}
	for i := 0; i < g.n; i++ {
		if negCol(g.col(w, i)) {
			out = append(out, i)
		}
	}
	return out
}

// word returns a reduced word by stripping the smallest right descent repeatedly.
func (g *En) word(w Elem) []int {
	w = g.clone(w)
	out := []int{}
	for {
		d := -1
		for i := 0; i < g.n; i++ {
			if negCol(g.col(w, i)) {
				d = i
				break
			}
		}
		if d < 0 {
			break
		}
		out = append(out, d)
		g.rmulInPlace(w, d)
	}
	for a, b := 0, len(out)-1; a < b; a, b = a+1, b-1 {
		out[a], out[b] = out[b], out[a]
	}
	return out
}

func (g *En) length(w Elem) int { return len(g.word(w)) }

func (g *En) inverse(w Elem) Elem {
	v := g.clone(g.e)
	word := g.word(w)
	for p := len(word) - 1; p >= 0; p-- {
		g.rmulInPlace(v, word[p])
	}
	return v
}

func (g *En) L(w Elem) []int { return g.R(g.inverse(w)) }

func (g *En) fromWord(word []int) Elem {
	w := g.clone(g.e)
	for _, i := range word {
		g.rmulInPlace(w, i)
	}
	return w
}

func (g *En) reflection(beta []int64) Elem {
	n := g.n
	out := make(Elem, n*n)
	for j := 0; j < n; j++ {
		m := g.pair(beta, j)
		for i := 0; i < n; i++ {
			var d int64
			if i == j {
				d = 1
			}
			out[j*n+i] = addChk(d, -mulChk(m, beta[i]))
		}
	}
	return out
}

func (g *En) rightTerminal(w Elem) bool {
	for _, s := range g.R(w) {
		for _, t := range g.adj[s] {
			sum := make([]int64, g.n)
			ws, wt := g.col(w, s), g.col(w, t)
			for k := range sum {
				sum[k] = addChk(ws[k], wt[k])
			}
			if !posCol(sum) {
				return false
			}
		}
	}
	return true
}

func (g *En) terminal(w Elem) bool {
	return g.rightTerminal(w) && g.rightTerminal(g.inverse(w))
}

// isFCWord is Stembridge's heap criterion for a reduced word in a simply laced
// group: fully commutative iff between two consecutive occurrences of a
// generator s there are at least two occurrences of neighbours of s.
func (g *En) isFCWord(word []int) bool {
	last := map[int]int{}
	for p, s := range word {
		if q, ok := last[s]; ok {
			cnt := 0
			for _, t := range word[q+1 : p] {
				if g.adjset[s][t] {
					cnt++
				}
			}
			if cnt < 2 {
				return false
			}
		}
		last[s] = p
	}
	return true
}

func elemKey(w Elem) string {
	var buf bytes.Buffer
	var b [8]byte
	for _, x := range w {
		u := uint64(x)
		for k := 0; k < 8; k++ {
			b[k] = byte(u >> (8 * k))
		}
		buf.Write(b[:])
	}
	return buf.String()
}

// covers returns all Bruhat covers of w (distinct elements) as reduced words:
// delete one letter of a reduced word, keep length l(w)-1.  Order of first
// occurrence.
func (g *En) covers(w Elem) []Elem {
	word := g.word(w)
	seen := map[string]bool{}
	out := []Elem{}
	for p := range word {
		sub := make([]int, 0, len(word)-1)
		sub = append(sub, word[:p]...)
		sub = append(sub, word[p+1:]...)
		v := g.fromWord(sub)
		key := elemKey(v)
		if !seen[key] && g.length(v) == len(word)-1 {
			seen[key] = true
			out = append(out, v)
		}
	}
	return out
}

// uniformData returns beta_r, gamma, delta of the manuscript's Section 4.1 in
// E_{4r+1}; a = 2r - 4.
func uniformData(r int) (n int, beta, gamma, delta []int64, a int64) {
	n = 4*r + 1
	beta = make([]int64, n)
	R := int64(r)
	beta[0], beta[1], beta[2], beta[4*r] = R-1, R, 2*R-1, R
	for j := 3; j < 4*r; j++ {
		beta[j] = 2*R - int64(j/2)
	}
	delta = make([]int64, n)
	gamma = make([]int64, n)
	for j, v := range []int64{2, 4, 6, 5, 4, 3, 2, 1} {
		delta[j] = v
	}
	delta[4*r] = 3
	for j, v := range []int64{1, 2, 3, 2, 2, 1, 1, 0} {
		gamma[j] = v
	}
	gamma[4*r] = 1
	return n, beta, gamma, delta, 2*R - 4
}

func betaRK(r, k int) []int64 {
	n, beta, gamma, delta, a := uniformData(r)
	K := int64(k)
	out := make([]int64, n)
	c1 := mulChk(a, K)
	c2 := addChk(mulChk(mulChk(a, K), K), mulChk(int64(r), K))
	for j := 0; j < n; j++ {
		out[j] = addChk(addChk(beta[j], -mulChk(c1, gamma[j])), mulChk(c2, delta[j]))
	}
	return out
}

func iR(r int) []int {
	set := map[int]bool{0: true, 4 * r: true}
	for j := 3; j < 4*r; j += 2 {
		set[j] = true
	}
	out := []int{}
	for j := range set {
		out = append(out, j)
	}
	sort.Ints(out)
	return out
}

func eqInts(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func eqElem(a, b Elem) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func contains(a []int, x int) bool {
	for _, y := range a {
		if y == x {
			return true
		}
	}
	return false
}

var failures = []string{}

func expect(cond bool, format string, args ...any) {
	if !cond {
		what := fmt.Sprintf(format, args...)
		failures = append(failures, what)
		fmt.Fprintln(os.Stderr, "EXPECTATION FAILED:", what)
	}
}

type bsRow struct {
	S                  int   `json:"s"`
	PairingM           int64 `json:"pairing_m"`
	BsLength           int   `json:"bs_length"`
	BsL                []int `json:"bs_L"`
	BsR                []int `json:"bs_R"`
	BsFullyCommutative bool  `json:"bs_fully_commutative"`
}

type pairRow struct {
	R                int     `json:"r"`
	K                int     `json:"k"`
	Rank             int     `json:"rank"`
	Beta             []int64 `json:"beta"`
	Ir               []int   `json:"I_r"`
	Length           int     `json:"length"`
	LengthFormula    int     `json:"length_formula_8r2_plus_3_plus_116k"`
	LEqualsREqualsIr bool    `json:"L_equals_R_equals_I_r"`
	Terminal         bool    `json:"terminal"`
	FullSupport      bool    `json:"full_support"`
	FullyCommutative bool    `json:"fully_commutative"`
	Involution       bool    `json:"involution"`
	BsElements       []bsRow `json:"bs_elements"`
}

type coverRow struct {
	R                   int `json:"r"`
	K                   int `json:"k"`
	Rank                int `json:"rank"`
	Length              int `json:"length"`
	DistinctCovers      int `json:"distinct_covers"`
	FullyCommutativeCov int `json:"fully_commutative_covers"`
	MinimumCoverSupport int `json:"minimum_cover_support"`
}

type report struct {
	Family             string     `json:"family"`
	Conventions        string     `json:"conventions"`
	CheckedPairs       [][2]int   `json:"checked_pairs"`
	LengthFormula      bool       `json:"length_formula_holds_for_all_checked_pairs"`
	BsChecks           string     `json:"bs_checks"`
	Pairs              []pairRow  `json:"pairs"`
	BruhatCovers       []coverRow `json:"bruhat_covers"`
	Status             string     `json:"status"`
	FailedExpectations []string   `json:"failed_expectations"`
}

func checkPair(r, k int) pairRow {
	n := 4*r + 1
	G := newEn(n)
	beta := betaRK(r, k)
	I := iR(r)
	expect(G.form(beta, beta) == 2, "r=%d k=%d: norm of beta is 2", r, k)
	allPos := true
	for _, c := range beta {
		if !(c > 0) {
			allPos = false
		}
	}
	expect(allPos, "r=%d k=%d: beta is positive", r, k)
	b := G.reflection(beta)
	word := G.word(b)
	predicted := 8*r*r + 3 + 116*k
	expect(len(word) == predicted, "r=%d k=%d: l(b) = 8r^2+3+116k", r, k)
	expect(len(word)%2 == 1, "r=%d k=%d: odd length", r, k)
	expect(eqInts(G.R(b), I) && eqInts(G.L(b), I), "r=%d k=%d: L(b) = R(b) = I_r", r, k)
	expect(len(I) == 2*r+1, "r=%d: |I_r| = 2r+1", r)
	indep := true
	for _, s := range I {
		for _, t := range G.adj[s] {
			if contains(I, t) {
				indep = false
			}
		}
	}
	expect(indep, "r=%d: I_r is independent", r)
	expect(G.terminal(b), "r=%d k=%d: terminal", r, k)
	support := map[int]bool{}
	for _, s := range word {
		support[s] = true
	}
	expect(len(support) == n, "r=%d k=%d: full support", r, k)
	expect(!G.isFCWord(word), "r=%d k=%d: b is not fully commutative", r, k)
	expect(eqElem(G.mul(b, b), G.e), "r=%d k=%d: involution", r, k)
	bsRows := []bsRow{}
	for _, s := range I {
		x := G.rmul(b, s)
		xw := G.word(x)
		m := G.pair(beta, s)
		sbeta := G.reflVec(beta, s)
		var rest []int
		for _, t := range I {
			if t != s {
				rest = append(rest, t)
			}
		}
		ok := len(xw) == len(word)-1 && eqInts(G.L(x), I) && eqInts(G.R(x), rest) && !G.isFCWord(xw)
		expect(ok, "r=%d k=%d s=%d: bs has L=I_r, R=I_r-{s}, length l(b)-1 and is not FC", r, k, s)
		expect(eqElem(G.inverse(x), G.lmul(s, b)), "r=%d k=%d s=%d: (bs)^-1 = sb", r, k, s)
		expect(m >= 1, "r=%d k=%d s=%d: (alpha_s, beta) >= 1 since s is a descent", r, k, s)
		nonneg := true
		var sum int64
		for _, c := range sbeta {
			if !(c >= 0) {
				nonneg = false
			}
			sum = addChk(sum, c)
		}
		expect(nonneg && sum > 1, "r=%d k=%d s=%d: s beta is a positive non-simple root", r, k, s)
		expect(G.pair(sbeta, s) == -m, "r=%d k=%d s=%d: (alpha_s, s beta) = -(alpha_s, beta)", r, k, s)
		expect(eqElem(G.lmul(s, x), G.reflection(sbeta)), "r=%d k=%d s=%d: s b s = r_(s beta)", r, k, s)
		bsRows = append(bsRows, bsRow{S: s, PairingM: m, BsLength: len(xw), BsL: G.L(x), BsR: G.R(x), BsFullyCommutative: G.isFCWord(xw)})
	}
	return pairRow{R: r, K: k, Rank: n, Beta: beta, Ir: I, Length: len(word), LengthFormula: predicted,
		LEqualsREqualsIr: true, Terminal: true, FullSupport: true, FullyCommutative: false, Involution: true, BsElements: bsRows}
}

func checkCovers(r, k int) coverRow {
	n := 4*r + 1
	G := newEn(n)
	b := G.reflection(betaRK(r, k))
	cov := G.covers(b)
	fc := 0
	minSupport := n
	for _, v := range cov {
		w := G.word(v)
		sup := map[int]bool{}
		for _, s := range w {
			sup[s] = true
		}
		if len(sup) < minSupport {
			minSupport = len(sup)
		}
		if G.isFCWord(w) {
			fc++
		}
	}
	expect(fc == 0, "r=%d k=%d: no Bruhat cover of b is fully commutative", r, k)
	return coverRow{R: r, K: k, Rank: n, Length: 8*r*r + 3 + 116*k, DistinctCovers: len(cov), FullyCommutativeCov: fc, MinimumCoverSupport: minSupport}
}

func main() {
	full := flag.Bool("full", false, "add r = 9, 10, 12 and further cover checks")
	flag.Parse()
	var pairs [][2]int
	for r := 3; r <= 8; r++ {
		for k := 0; k < 4; k++ {
			pairs = append(pairs, [2]int{r, k})
		}
	}
	coverPairs := [][2]int{{3, 0}, {3, 1}, {4, 0}}
	if *full {
		pairs = append(pairs, [2]int{9, 0}, [2]int{9, 3}, [2]int{10, 0}, [2]int{10, 2}, [2]int{12, 0}, [2]int{12, 1})
		coverPairs = append(coverPairs, [2]int{3, 2}, [2]int{4, 1}, [2]int{5, 0}, [2]int{6, 0}, [2]int{7, 0})
	}
	rows := []pairRow{}
	for _, p := range pairs {
		rows = append(rows, checkPair(p[0], p[1]))
	}
	covers := []coverRow{}
	for _, p := range coverPairs {
		covers = append(covers, checkCovers(p[0], p[1]))
	}
	expected := map[[2]int]int{{3, 0}: 74, {3, 1}: 162, {3, 2}: 230, {4, 0}: 130, {4, 1}: 218, {5, 0}: 202, {6, 0}: 290, {7, 0}: 394}
	for _, c := range covers {
		key := [2]int{c.R, c.K}
		if want, ok := expected[key]; ok {
			expect(c.DistinctCovers == want, "r=%d k=%d: %d distinct covers", key[0], key[1], want)
		}
	}
	lengthOK := true
	for _, r := range rows {
		if r.Length != r.LengthFormula {
			lengthOK = false
		}
	}
	status := "passed"
	if len(failures) > 0 {
		status = "FAILED"
	}
	rep := report{
		Family:        "b_{r,k} = r_{beta_{r,k}} in W(E_{4r+1}), beta_{r,k} = beta_r - a k gamma + (a k^2 + r k) delta, a = 2r-4",
		Conventions:   "E_n: chain 0-1-...-(n-2), node n-1 attached to node 2; words act on the right; columns of w are w(alpha_j)",
		CheckedPairs:  pairs,
		LengthFormula: lengthOK,
		BsChecks:      "for every s in I_r: L(bs) = I_r, R(bs) = I_r - {s}, l(bs) = l(b) - 1, bs not FC, sbs = r_{s beta} with s beta positive non-simple and (alpha_s, s beta) = -(alpha_s, beta) < 0",
		Pairs:         rows, BruhatCovers: covers, Status: status, FailedExpectations: failures,
	}
	outPath := filepath.Join("results", "uniform-family-certificate.json")
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", " ")
	if err := enc.Encode(rep); err != nil {
		fmt.Fprintln(os.Stderr, "encode:", err)
		os.Exit(2)
	}
	if err := os.MkdirAll("results", 0o755); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}
	if err := os.WriteFile(outPath, buf.Bytes(), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}
	// One-line summary on stdout, in the format of the original.
	var sb strings.Builder
	fmt.Fprintf(&sb, `{"status": %q, "pairs_checked": %d, "lengths": {`, status, len(rows))
	for i, r := range rows {
		if i > 0 {
			sb.WriteString(", ")
		}
		fmt.Fprintf(&sb, `"(%d,%d)": %d`, r.R, r.K, r.Length)
	}
	sb.WriteString(`}, "covers": {`)
	for i, c := range covers {
		if i > 0 {
			sb.WriteString(", ")
		}
		fmt.Fprintf(&sb, `"(%d,%d)": [%d, %d]`, c.R, c.K, c.DistinctCovers, c.FullyCommutativeCov)
	}
	sb.WriteString("}}")
	fmt.Println(sb.String())
	if len(failures) > 0 {
		os.Exit(1)
	}
}
