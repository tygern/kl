// Command uniform-construction runs the finite counterchecks of the
// rank-uniform E_(4r+1) construction (r >= 3). The all-r proof is in
// research/en_uniform/construction.txt; the finite samples here do not replace
// it. For each r in 3..max-r and k in {0,1,2,10} it verifies, with exact
// integer arithmetic on simple-root coordinates, the signed-coordinate seed
// certificate, the translation identities T(e_i), T(b), T(v), T(u), the
// sign pattern of the pairing coefficients, and for every k that the root
// beta_(r,k) is a norm-2 positive root whose reflection has right descent set
// I and is right terminal (hence, being an involution, also left terminal).
// No group-ball, Bruhat-interval or FC-catalogue enumeration is used.
//
// Writes research/en_uniform/construction.json (relative to the current
// working directory, overridable with -out) and prints the same two summary
// lines as the original.
//
// Ports research/en_uniform/construction.py (with the E_n Coxeter helper
// research/en_families/targeted.py: only the adjacency, right-multiplication
// on matrix-of-columns elements and descent set are needed, ported inline).
//
// Imports: only the Go standard library; no package of this module.
// Run: uniform-construction -max-r 30
package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
)

const defaultOut = "research/en_uniform/construction.json"

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "uniform-construction: assertion failed: "+format+"\n", args...)
	os.Exit(1)
}

func must(cond bool, what string) {
	if !cond {
		fail("%s", what)
	}
}

// Checked int64 arithmetic: any overflow aborts.
func addc(a, b int64) int64 {
	c := a + b
	if (b > 0 && c < a) || (b < 0 && c > a) {
		panic("int64 overflow in addition")
	}
	return c
}

func mulc(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	c := a * b
	if c/b != a {
		panic("int64 overflow in multiplication")
	}
	return c
}

type vec = []int64

// cox is the simply-laced Coxeter system E_n (nodes 0..n-1; path 0-1-...-(n-2)
// plus the edge 2-(n-1)).
type cox struct {
	n   int
	adj [][]int
	e   []vec
}

func en(n int) *cox {
	g := &cox{n: n, adj: make([][]int, n)}
	add := func(s, t int) {
		g.adj[s] = append(g.adj[s], t)
		g.adj[t] = append(g.adj[t], s)
	}
	for i := 0; i < n-2; i++ {
		add(i, i+1)
	}
	add(2, n-1)
	g.e = make([]vec, n)
	for i := range g.e {
		g.e[i] = make(vec, n)
		g.e[i][i] = 1
	}
	return g
}

// right is right multiplication of the matrix-of-columns element w by the
// simple reflection s (rows of w are the images of the simple roots).
func (g *cox) right(w []vec, s int) []vec {
	out := make([]vec, len(w))
	copy(out, w)
	neg := make(vec, g.n)
	for i, a := range w[s] {
		neg[i] = -a
	}
	out[s] = neg
	for _, t := range g.adj[s] {
		row := make(vec, g.n)
		for i := range row {
			row[i] = addc(w[s][i], w[t][i])
		}
		out[t] = row
	}
	return out
}

func (g *cox) desc(w []vec) []int {
	out := []int{}
	for i, row := range w {
		all := true
		for _, x := range row {
			if x > 0 {
				all = false
				break
			}
		}
		if all {
			out = append(out, i)
		}
	}
	return out
}

func (g *cox) pair(x, y vec) int64 {
	var sum int64
	for i := 0; i < g.n; i++ {
		if x[i] == 0 {
			continue
		}
		inner := mulc(2, y[i])
		for _, j := range g.adj[i] {
			inner = addc(inner, -y[j])
		}
		sum = addc(sum, mulc(x[i], inner))
	}
	return sum
}

func (g *cox) reflect(root, x vec) vec {
	m := g.pair(root, x)
	out := make(vec, g.n)
	for i := range out {
		out[i] = addc(x[i], -mulc(m, root[i]))
	}
	return out
}

func sum(x vec) int64 {
	var s int64
	for _, v := range x {
		s = addc(s, v)
	}
	return s
}

func minOf(x vec) int64 {
	m := x[0]
	for _, v := range x {
		if v < m {
			m = v
		}
	}
	return m
}

func maxOf(x vec) int64 {
	m := x[0]
	for _, v := range x {
		if v > m {
			m = v
		}
	}
	return m
}

func equal(a, b vec) bool {
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

func contains(s []int, x int) bool {
	for _, v := range s {
		if v == x {
			return true
		}
	}
	return false
}

func equalInts(a, b []int) bool {
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

func seed(r int) vec {
	out := vec{int64(r - 1), int64(r), int64(2*r - 1)}
	for j := 3; j < 4*r; j++ {
		out = append(out, int64((4*r+1-j)/2)) // 4r+1-j > 0: floor == truncation
	}
	return append(out, int64(r))
}

type check struct {
	K                 int64 `json:"k"`
	Height            int64 `json:"height"`
	LargestCoordinate int64 `json:"largest_coordinate"`
	Norm              int   `json:"norm"`
	Descents          []int `json:"descents"`
	Terminal          bool  `json:"terminal"`
}

type row struct {
	R                       int     `json:"r"`
	Rank                    int     `json:"rank"`
	Seed                    vec     `json:"seed"`
	I                       []int   `json:"I"`
	Delta                   vec     `json:"delta"`
	Gamma                   vec     `json:"gamma"`
	A                       int64   `json:"a"`
	SeedS0DCoordinatesTwice vec     `json:"seed_s0_D_coordinates_twice"`
	EvenSignChange          []int   `json:"even_sign_change_coordinates"`
	PairingCoefficients     []vec   `json:"pairing_coefficients"`
	Checks                  []check `json:"checks"`
}

type output struct {
	Scope            string `json:"scope"`
	TheoremRMinimum  int    `json:"theorem_r_minimum"`
	CheckedRMaximum  int    `json:"checked_r_maximum"`
	ParameterFormula string `json:"parameter_formula"`
	Rows             []row  `json:"rows"`
}

func verify(r int, ks []int64) row {
	n := 4*r + 1
	N := 4 * r
	g := en(n)
	b := seed(r)
	a := int64(2*r - 4)
	rr := int64(r)
	I := []int{0}
	for j := 3; j < 4*r; j += 2 {
		I = append(I, j)
	}
	I = append(I, 4*r)
	delta := make(vec, n)
	copy(delta, vec{2, 4, 6, 5, 4, 3, 2, 1})
	delta[n-1] = 3
	gamma := make(vec, n)
	copy(gamma, vec{1, 2, 3, 2, 2, 1, 1, 0})
	gamma[n-1] = 1

	// Twice the finite-D component in alpha0=h+lambda coordinates.
	dcoordsTwice := func(c vec) vec {
		out := vec{-c[0] + 2*c[1] - 2*c[N], c[0] - 2*c[1] + 2*c[2] - 2*c[N]}
		for j := 3; j < N; j++ {
			out = append(out, c[0]-2*c[j-1]+2*c[j])
		}
		return append(out, c[0]-2*c[N-1])
	}

	afterS0 := g.reflect(g.e[0], b)
	must(afterS0[0] == 1 && equal(afterS0[1:], b[1:]), "reflection of seed in s0")
	lambdaTwice := make(vec, 0, N)
	lambdaTwice = append(lambdaTwice, -1)
	for i := 1; i < N; i++ {
		lambdaTwice = append(lambdaTwice, 1)
	}
	signs := make(vec, 0, N)
	negs := 0
	for j := 1; j <= N; j++ {
		if j%2 == 0 {
			signs = append(signs, -1)
			negs++
		} else {
			signs = append(signs, 1)
		}
	}
	must(negs == 2*r, "even signed permutation in W(D_N)") // even signed permutation, in W(D_N)
	dc := dcoordsTwice(afterS0)
	want := make(vec, N)
	for i := range want {
		want[i] = signs[i] * lambdaTwice[i]
	}
	must(equal(dc, want), "seed D-coordinates")
	must(g.pair(b, b) == 2 && minOf(b) > 0, "seed norm and positivity")
	must(g.pair(gamma, gamma) == 2 && g.pair(delta, delta) == 0 && g.pair(gamma, delta) == 0, "gamma, delta norms")
	must(g.pair(b, gamma) == -rr && g.pair(b, delta) == -a, "pairings of seed with gamma, delta")
	gd := make(vec, n)
	for i := range gd {
		gd[i] = gamma[i] + delta[i]
	}
	gammaSteps := []int{2, 1, 0, 4, 3, 2, 1, 6, 5, 4, 3, 2}
	gdSteps := append(append([]int{}, gammaSteps...),
		N, 2, 1, 0, 3, 2, 1, 4, 3, 2, 5, 4, 3, 6, 5, 4, 7, 6, 5, N, 2, 1, 0, 3, 2, 1, 4, 3, 2)
	type pr struct {
		root  vec
		steps []int
	}
	for _, p := range []pr{{gamma, gammaSteps}, {gd, gdSteps}} {
		root := p.root
		for _, s := range p.steps {
			oldHeight := sum(root)
			root = g.reflect(g.e[s], root)
			must(minOf(root) >= 0 && sum(root) < oldHeight, "height descent of gamma, gamma+delta")
		}
		must(equal(root, g.e[N]), "gamma, gamma+delta reduce to alpha_N")
	}

	T := func(x vec) vec { return g.reflect(gamma, g.reflect(gd, x)) }

	for _, e := range g.e {
		pe := g.pair(e, delta)
		pg := g.pair(e, gamma)
		got := T(e)
		for i := 0; i < n; i++ {
			exp := addc(addc(e[i], mulc(pe, gamma[i])), -mulc(addc(pg, pe), delta[i]))
			must(got[i] == exp, "T on simple roots")
		}
	}
	v := make(vec, n)
	u := make(vec, n)
	for i := range v {
		v[i] = addc(mulc(a+rr, delta[i]), -mulc(a, gamma[i]))
		u[i] = mulc(2*a, delta[i])
	}
	anyV := false
	for _, x := range v {
		if x != 0 {
			anyV = true
		}
	}
	must(minOf(v) >= 0 && minOf(u) >= 0 && anyV, "v, u nonnegative and v nonzero")
	plus := func(x, y vec) vec {
		out := make(vec, n)
		for i := range out {
			out[i] = addc(x[i], y[i])
		}
		return out
	}
	must(equal(T(b), plus(b, v)), "T(b)=b+v")
	must(equal(T(v), plus(v, u)) && equal(T(u), u), "T(v)=v+u, T(u)=u")
	coeffs := make([]vec, 0, 3)
	for _, z := range []vec{b, v, u} {
		m := make(vec, n)
		for i, e := range g.e {
			m[i] = g.pair(e, z)
		}
		coeffs = append(coeffs, m)
	}
	inI := make([]bool, n)
	for _, i := range I {
		inI[i] = true
	}
	for degree, m := range coeffs {
		lo := int64(0)
		if degree == 0 {
			lo = 1
		}
		for _, i := range I {
			must(m[i] >= lo, "pairing coefficient lower bound on I")
		}
		for i := 0; i < n; i++ {
			if !inI[i] {
				must(m[i] <= 0, "pairing coefficient sign off I")
			}
		}
		for _, s := range I {
			for _, t := range g.adj[s] {
				must(m[s]+m[t] <= 0, "adjacent pairing sums")
			}
		}
	}
	matching := [][2]int{{0, 1}, {2, N}}
	for j := 3; j < N-1; j += 2 {
		matching = append(matching, [2]int{j, j + 1})
	}
	seen := map[int]bool{}
	for _, e := range matching {
		seen[e[0]] = true
		seen[e[1]] = true
	}
	must(len(matching) == 2*r && len(seen) == 4*r, "matching")
	must(len(I) == 2*r+1, "size of I")
	for _, s := range I {
		for _, t := range g.adj[s] {
			must(!inI[t], "I independent")
		}
	}
	checks := []check{}
	for _, k := range ks {
		root := make(vec, n)
		for i := range root {
			root[i] = addc(addc(b[i], -mulc(mulc(a, k), gamma[i])), mulc(addc(mulc(mulc(a, k), k), mulc(rr, k)), delta[i]))
		}
		expected := vec{rr - 2, 2 - rr, -1 - a*k}
		for j := 3; j < N; j++ {
			if (j+1)%2 == 0 {
				expected = append(expected, 1)
			} else {
				expected = append(expected, -1)
			}
		}
		expected = append(expected, 1+a*k)
		for j := 3; j < 8; j++ {
			if (j+1)%2 == 0 {
				expected[j] = 1 + a*k
			} else {
				expected[j] = -(1 + a*k)
			}
		}
		expected[8] = -1 - a*k*k - rr*k
		must(minOf(root) > 0 && g.pair(root, root) == 2, "root positive of norm 2")
		got := make(vec, n)
		for i, e := range g.e {
			got[i] = g.pair(e, root)
		}
		must(equal(got, expected), "pairing of root with simple roots")
		w := make([]vec, n)
		for i, e := range g.e {
			w[i] = g.reflect(root, e)
		}
		must(equalInts(g.desc(w), I), "descent set of reflection")
		// Right terminality and involution imply left terminality without
		// reconstructing increasingly long reduced words.
		for _, s := range I {
			ds := g.desc(g.right(w, s))
			for _, t := range g.adj[s] {
				must(!contains(ds, t), "right terminality")
			}
		}
		for _, e := range g.e {
			must(equal(g.reflect(root, g.reflect(root, e)), e), "reflection is an involution")
		}
		checks = append(checks, check{K: k, Height: sum(root), LargestCoordinate: maxOf(root), Norm: 2, Descents: I, Terminal: true})
	}
	even := []int{}
	for j := 2; j <= N; j += 2 {
		even = append(even, j)
	}
	return row{R: r, Rank: n, Seed: b, I: I, Delta: delta, Gamma: gamma, A: a,
		SeedS0DCoordinatesTwice: dcoordsTwice(afterS0), EvenSignChange: even,
		PairingCoefficients: coeffs, Checks: checks}
}

func main() {
	maxR := flag.Int("max-r", 30, "largest r to check (at least 3)")
	out := flag.String("out", defaultOut, "output JSON path (relative to the working directory)")
	flag.Parse()
	if flag.NArg() != 0 || *maxR < 3 {
		fmt.Fprintln(os.Stderr, "usage: uniform-construction [-max-r N] [-out PATH]   (N >= 3; no positional arguments)")
		os.Exit(2)
	}
	rows := []row{}
	for r := 3; r <= *maxR; r++ {
		rows = append(rows, verify(r, []int64{0, 1, 2, 10}))
	}
	o := output{
		Scope:            "All-r proof is construction.txt; these are finite counterchecks",
		TheoremRMinimum:  3,
		CheckedRMaximum:  *maxR,
		ParameterFormula: "beta_(r,k)=beta_r-(2r-4)k*gamma+((2r-4)k^2+r*k)*delta",
		Rows:             rows,
	}
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(o); err != nil {
		fail("encoding: %v", err)
	}
	if dir := filepath.Dir(*out); dir != "." {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			fail("creating directory: %v", err)
		}
	}
	if err := os.WriteFile(*out, buf.Bytes(), 0o644); err != nil {
		fail("writing %s: %v", *out, err)
	}
	fmt.Printf("Verified r=3,...,%d; ranks13,...,%d; k=0,1,2,10.\n", *maxR, 4**maxR+1)
	fmt.Println("Checked finite-D signed-coordinate seed certificate and symbolic translation identities.")
}
