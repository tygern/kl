// Command uniform-verify is the independent exact audit of the rank-uniform
// E_(4r+1) family. It checks, in the exact polynomial ring Q[r,k,q], the local
// identities of the all-rank proof (math/big.Rat coefficients), and then, for
// each r=3..40 and k in {0,1,2,10,10^6}, exact integer matrix and root
// checks. It writes research/en_uniform/referee-certificate.json (relative to
// the current working directory), including the SHA-256 of
// results/exceptional-leading.tex, and prints the same three summary lines as
// the original.
//
// Ports research/en_uniform/referee_verify.py (the Python Fraction class
// becomes math/big.Rat; Python integers become int64 with checked arithmetic
// (roots reach 6e13) and math/big.Int for the reflection columns, which reach
// about 1e22).
//
// Imports: only the Go standard library. It imports no package of this
// module: no construction, targeted, or FC enumeration code. Run: uniform-verify
package main

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"math/big"
	"os"
	"reflect"
)

const (
	sourcePath = "results/exceptional-leading.tex"
	outputPath = "research/en_uniform/referee-certificate.json"
)

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "uniform-verify: assertion failed: "+format+"\n", args...)
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
	if c/b != a || (a == -1 && b == math.MinInt64) || (b == -1 && a == math.MinInt64) {
		panic("int64 overflow in multiplication")
	}
	return c
}

type symbolicOut struct {
	CoefficientRing string `json:"coefficient_ring"`
	AllPass         bool   `json:"all_assertions_pass"`
	Edge78Sum       string `json:"edge_7_8_sum"`
	Edge89Sum       string `json:"edge_8_9_sum"`
}

func symbolic() symbolicOut {
	r, k, q := variable(0), variable(1), variable(2)
	c := cst
	half := frac(1, 2)
	a := c(2).mul(r).sub(c(4))
	// Seed norm after cancelling paired terms on the long arm.
	must(r.sub(c(1)).mul(r.sub(c(2))).add(r.mul(c(2).sub(r))).add(r).eq(c(2)), "seed norm")
	// In the independently chosen D-coordinate realization, s_0 beta has
	// alpha_0 coefficient one and alternating finite coordinates -1/2,+1/2.
	must(r.sub(c(1)).sub(r.sub(c(2))).eq(c(1)), "alpha_0 coefficient")
	must(r.sub(c(1)).neg().mul(half).add(r.sub(c(2)).mul(half)).eq(frac(-1, 2)), "coordinate -1/2")
	must(c(3).sub(r).mul(half).add(r.sub(c(2)).mul(half)).eq(frac(1, 2)), "coordinate +1/2")
	// Generic interior tail pairing, at even 2q and odd 2q+1.
	even := c(2).mul(r).sub(q)
	odd := c(2).mul(r).sub(q)
	must(c(2).mul(even).sub(c(2).mul(r).sub(q).add(c(1))).sub(odd).eq(c(-1)), "tail pairing even")
	must(c(2).mul(odd).sub(even).sub(c(2).mul(r).sub(q).sub(c(1))).eq(c(1)), "tail pairing odd")
	// The only nonconstant local section is supported in nodes 0..8,branch.
	delta := [10]int64{2, 4, 6, 5, 4, 3, 2, 1, 0, 0}
	gamma := [10]int64{1, 2, 3, 2, 2, 1, 1, 0, 0, 0}
	seed := []poly{r.sub(c(1)), r}
	for j := 2; j < 10; j++ {
		seed = append(seed, c(2).mul(r).sub(c(int64(j/2))))
	}
	ak := a.mul(k)
	kk := a.mul(k).mul(k).add(r.mul(k)) // a*k*k + r*k
	beta := make([]poly, 10)
	for j := 0; j < 10; j++ {
		beta[j] = seed[j].sub(ak.mul(c(gamma[j]))).add(kk.mul(c(delta[j])))
	}
	branch := r.sub(ak).add(kk.mul(c(3)))
	pairings := []poly{
		c(2).mul(beta[0]).sub(beta[1]),
		c(2).mul(beta[1]).sub(beta[0]).sub(beta[2]),
		c(2).mul(beta[2]).sub(beta[1]).sub(beta[3]).sub(branch),
	}
	for j := 3; j < 9; j++ {
		pairings = append(pairings, c(2).mul(beta[j]).sub(beta[j-1]).sub(beta[j+1]))
	}
	one := c(1)
	expected := []poly{r.sub(c(2)), c(2).sub(r), c(-1).sub(ak)}
	for j := 3; j < 8; j++ {
		sign := int64(1)
		if (j+1)%2 != 0 {
			sign = -1
		}
		expected = append(expected, c(sign).mul(one.add(ak)))
	}
	expected = append(expected, c(-1).sub(a.mul(k).mul(k)).sub(r.mul(k)))
	must(len(pairings) == len(expected), "pairings length")
	for j := range pairings {
		must(pairings[j].eq(expected[j]), fmt.Sprintf("pairing %d", j))
	}
	must(c(2).mul(branch).sub(beta[2]).eq(one.add(ak)), "branch pairing")
	m0 := []poly{r.sub(c(2)), c(2).sub(r), c(-1), c(1), c(-1), c(1), c(-1), c(1)}
	sg, sd := c(0), c(0)
	for j := 0; j < 8; j++ {
		sg = sg.add(m0[j].mul(c(gamma[j])))
		sd = sd.add(m0[j].mul(c(delta[j])))
	}
	must(sg.add(c(1)).eq(r.neg()), "gamma pairing -r")
	must(sd.add(c(3)).eq(a.neg()), "delta pairing -a")
	// Exceptional edge inequalities have nonpositive binomial coefficients.
	kkm1 := a.mul(k).mul(k.sub(c(1))) // a*k*(k-1)
	must(expected[7].add(expected[8]).eq(r.neg().mul(k).sub(kkm1)), "edge 7-8")
	must(expected[8].add(c(1)).eq(a.add(r).neg().mul(k).sub(kkm1)), "edge 8-9")
	must(beta[0].eq(r.sub(c(1)).add(c(4).mul(k)).add(c(4).mul(r).sub(c(8)).mul(k).mul(k))), "beta0")
	return symbolicOut{
		CoefficientRing: "Q[r,k,q]",
		AllPass:         true,
		Edge78Sum:       "-r*k-(2*r-4)*k*(k-1)",
		Edge89Sum:       "-(3*r-4)*k-(2*r-4)*k*(k-1)",
	}
}

type witness struct {
	LoweringWord []any `json:"lowering_word"`
	EndsAt       any   `json:"ends_at"`
}

type sample struct {
	K           int64 `json:"k"`
	Height      int64 `json:"height"`
	Coordinate0 int64 `json:"coordinate_0"`
}

type rowOut struct {
	R            int      `json:"r"`
	Rank         int      `json:"rank"`
	Alpha        int      `json:"alpha"`
	SampleChecks []sample `json:"sample_checks"`
}

type vec = []int64

func vecEq(x, y vec) bool {
	if len(x) != len(y) {
		return false
	}
	for i := range x {
		if x[i] != y[i] {
			return false
		}
	}
	return true
}

func vecSum(x vec) int64 {
	var s int64
	for _, v := range x {
		s = addc(s, v)
	}
	return s
}

func vecMin(x vec) int64 {
	m := x[0]
	for _, v := range x {
		if v < m {
			m = v
		}
	}
	return m
}

func vecMax(x vec) int64 {
	m := x[0]
	for _, v := range x {
		if v > m {
			m = v
		}
	}
	return m
}

func checkRank(r int) (rowOut, []witness) {
	n := 4*r + 1
	a := int64(2*r - 4)
	ri := int64(r)
	type edge struct{ s, t int }
	var edges []edge
	for j := 0; j < n-2; j++ {
		edges = append(edges, edge{j, j + 1})
	}
	edges = append(edges, edge{2, n - 1})
	adj := make([][]int, n)
	for _, e := range edges {
		adj[e.s] = append(adj[e.s], e.t)
		adj[e.t] = append(adj[e.t], e.s)
	}
	basis := make([]vec, n)
	for j := 0; j < n; j++ {
		basis[j] = make(vec, n)
		basis[j][j] = 1
	}
	beta0 := vec{ri - 1, ri}
	for j := 2; j < n-1; j++ {
		beta0 = append(beta0, int64(2*r-j/2))
	}
	beta0 = append(beta0, ri)
	delta := vec{2, 4, 6, 5, 4, 3, 2, 1}
	gamma := vec{1, 2, 3, 2, 2, 1, 1, 0}
	for i := 0; i < n-9; i++ {
		delta = append(delta, 0)
		gamma = append(gamma, 0)
	}
	delta = append(delta, 3)
	gamma = append(gamma, 1)
	I := map[int]bool{0: true, n - 1: true}
	for j := 3; j < n-1; j += 2 {
		I[j] = true
	}
	inI := func(j int) bool { return I[j] }

	cartan := func(x vec) vec {
		out := make(vec, n)
		for j := 0; j < n; j++ {
			s := mulc(2, x[j])
			for _, t := range adj[j] {
				s = addc(s, -x[t])
			}
			out[j] = s
		}
		return out
	}
	pair := func(x, y vec) int64 {
		// Individual products reach about 1e22 for the k=10^6 roots, so
		// accumulate exactly in big.Int.
		cy := cartan(y)
		s := new(big.Int)
		for i := range x {
			s.Add(s, new(big.Int).Mul(big.NewInt(x[i]), big.NewInt(cy[i])))
		}
		if !s.IsInt64() {
			panic("pairing does not fit int64")
		}
		return s.Int64()
	}
	add := func(x, y vec, scale int64) vec {
		out := make(vec, n)
		for i := range x {
			out[i] = addc(x[i], mulc(scale, y[i]))
		}
		return out
	}
	refl := func(root, x vec) vec {
		return add(x, root, -pair(root, x))
	}
	gd := add(gamma, delta, 1)
	T := func(x vec) vec { return refl(gamma, refl(gd, x)) }
	N := func(x vec) vec { return add(T(x), x, -1) }
	zero := make(vec, n)
	basisIndex := func(x vec) int {
		for i, b := range basis {
			if vecEq(x, b) {
				return i
			}
		}
		return -1
	}

	// Uniform positive-real-root witness for gamma and gamma+delta, supported
	// in the same affine parabolic. The lowering path certifies real-root
	// membership by reduction to a simple root. Norm two alone is insufficient.
	var witnesses []witness
	for _, root := range []vec{gamma, gd} {
		x := root
		word := []any{}
		for vecSum(x) > 1 {
			p := cartan(x)
			s := -1
			for j := 0; j < n; j++ {
				if p[j] > 0 {
					s = j
					break
				}
			}
			must(s >= 0, "no positive Cartan coordinate in lowering")
			y := refl(basis[s], x)
			must(vecMin(y) >= 0 && vecSum(y) < vecSum(x), "lowering step")
			must(s < 8 || s == n-1, "lowering node in affine parabolic")
			if s == n-1 {
				word = append(word, "branch")
			} else {
				word = append(word, s)
			}
			x = y
		}
		must(basisIndex(x) >= 0 && pair(root, root) == 2, "witness ends at simple root of norm 2")
		var ends any
		if vecEq(x, basis[n-1]) {
			ends = "branch"
		} else {
			ends = basisIndex(x)
		}
		witnesses = append(witnesses, witness{LoweringWord: word, EndsAt: ends})
	}
	must(pair(delta, delta) == 0 && pair(gamma, delta) == 0, "delta null and orthogonal to gamma")
	must(pair(beta0, gamma) == -ri && pair(beta0, delta) == -a, "beta0 pairings")
	for _, e := range basis {
		want := add(add(e, gamma, pair(e, delta)), delta, -pair(e, gamma)-pair(e, delta))
		must(vecEq(T(e), want), "T on basis")
		must(vecEq(N(N(N(e))), zero), "N cubed zero")
	}
	v := make(vec, n)
	u := make(vec, n)
	for i := 0; i < n; i++ {
		v[i] = addc(mulc(a+ri, delta[i]), -mulc(a, gamma[i]))
		u[i] = mulc(2*a, delta[i])
	}
	must(vecEq(N(beta0), v) && vecEq(N(v), u) && vecEq(N(u), zero), "N chain")
	must(vecMin(beta0) > 0 && vecMin(v) >= 0 && vecMin(u) >= 0, "positivity of beta0, v, u")
	// Independent D-coordinate realization, doubled to avoid fractions.
	after0 := refl(basis[0], beta0)
	finite := make(vec, n-1)
	for i := range finite {
		finite[i] = -after0[0]
	}
	finite[0] += 2*after0[1] - 2*after0[n-1]
	finite[1] += 2*after0[1] + 2*after0[n-1]
	for j := 2; j < n-1; j++ {
		finite[j] += 2 * after0[j]
		finite[j-1] -= 2 * after0[j]
	}
	must(after0[0] == 1, "s_0 beta alpha_0 coefficient")
	wantFinite := make(vec, 0, n-1)
	for i := 1; i < n; i++ {
		if i%2 == 1 {
			wantFinite = append(wantFinite, -1)
		} else {
			wantFinite = append(wantFinite, 1)
		}
	}
	must(vecEq(finite, wantFinite), "alternating finite coordinates")
	npos := 0
	for _, cval := range finite {
		if cval > 0 {
			npos++
		}
	}
	must(npos == 2*r, "allowed even sign changes")
	matching := [][2]int{{0, 1}, {2, n - 1}}
	for j := 3; j < n-2; j += 2 {
		matching = append(matching, [2]int{j, j + 1})
	}
	covered := map[int]bool{}
	for _, e := range matching {
		covered[e[0]] = true
		covered[e[1]] = true
	}
	must(len(matching) == 2*r && len(covered) == 4*r, "matching")
	must(len(I) == n-len(matching) && n-len(matching) == 2*r+1, "size of I")
	for _, e := range edges {
		must(!(inI(e.s) && inI(e.t)), "I independent")
	}
	m0, m1, m2 := cartan(beta0), cartan(v), cartan(u)
	for j := 0; j < n; j++ {
		if inI(j) {
			must(m0[j] > 0 && m1[j] >= 0 && m2[j] >= 0, "signs on I")
		} else {
			must(m0[j] < 0 && m1[j] <= 0 && m2[j] <= 0, "signs off I")
		}
	}
	for _, e := range edges {
		if inI(e.s) || inI(e.t) {
			for _, m := range []vec{m0, m1, m2} {
				must(m[e.s]+m[e.t] <= 0, "edge sums nonpositive")
			}
		}
	}
	var rows []sample
	for _, k := range []int64{0, 1, 2, 10, 1000000} {
		root := make(vec, n)
		for i := 0; i < n; i++ {
			t := addc(beta0[i], -mulc(mulc(a, k), gamma[i]))
			coef := addc(mulc(mulc(a, k), k), mulc(ri, k))
			root[i] = addc(t, mulc(coef, delta[i]))
		}
		// For k >= 0, k*(k-1) >= 0 and floor and truncating division return the same integer.
		tri := mulc(k, k-1) / 2
		for i := 0; i < n; i++ {
			alt := addc(addc(beta0[i], mulc(k, v[i])), mulc(tri, u[i]))
			must(root[i] == alt, "root equals quadratic family in v,u")
		}
		must(vecMin(root) > 0 && pair(root, root) == 2, "root positive of norm 2")
		m := cartan(root)
		for j := 0; j < n; j++ {
			must((m[j] > 0) == inI(j), "positive Cartan coordinates are exactly I")
		}
		// Reflection columns are e_j - pair(root,e_j)*root. Their entries reach
		// about 1e22 at k=10^6, beyond int64, so they are exact big integers.
		columns := make([][]*big.Int, n)
		for j, e := range basis {
			sc := big.NewInt(pair(root, e))
			col := make([]*big.Int, n)
			for i := 0; i < n; i++ {
				col[i] = new(big.Int).Sub(big.NewInt(e[i]), new(big.Int).Mul(sc, big.NewInt(root[i])))
			}
			columns[j] = col
		}
		for j, col := range columns {
			mx := col[0]
			for _, c := range col {
				if c.Cmp(mx) > 0 {
					mx = c
				}
			}
			must((mx.Sign() <= 0) == inI(j), "nonpositive columns are exactly I")
		}
		for _, e := range edges {
			if inI(e.s) || inI(e.t) {
				for i := 0; i < n; i++ {
					must(new(big.Int).Add(columns[e.s][i], columns[e.t][i]).Sign() >= 0, "edge column sums nonnegative")
				}
			}
		}
		for j := 0; j < n; j++ {
			must(columns[0][j].Cmp(big.NewInt(basis[0][j])) != 0, "column 0 differs in every coordinate")
		}
		must(root[0] == addc(addc(ri-1, mulc(4, k)), mulc(mulc(int64(4*r-8), k), k)), "coordinate 0 formula")
		rows = append(rows, sample{K: k, Height: vecSum(root), Coordinate0: root[0]})
	}
	return rowOut{R: r, Rank: n, Alpha: len(I), SampleChecks: rows}, witnesses
}

type certificate struct {
	Symbolic     symbolicOut `json:"symbolic"`
	Scope        string      `json:"scope"`
	Independence string      `json:"independence"`
	Rows         []rowOut    `json:"rows"`
	Witnesses    []witness   `json:"uniform_affine_root_witnesses"`
	SourceSHA256 string      `json:"audited_source_sha256"`
}

func main() {
	data := certificate{
		Symbolic:     symbolic(),
		Scope:        "All-rank proof in referee.txt; matrix rows are finite counterchecks.",
		Independence: "No construction, targeted, or FC enumeration code is imported.",
		Rows:         []rowOut{},
	}
	for r := 3; r <= 40; r++ {
		row, witnesses := checkRank(r)
		if r == 3 {
			data.Witnesses = witnesses
		} else {
			must(reflect.DeepEqual(witnesses, data.Witnesses), "witnesses identical across ranks")
		}
		data.Rows = append(data.Rows, row)
	}
	src, err := os.ReadFile(sourcePath)
	if err != nil {
		fail("cannot read %s: %v", sourcePath, err)
	}
	sum := sha256.Sum256(src)
	data.SourceSHA256 = hex.EncodeToString(sum[:])
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(data); err != nil {
		fail("encode: %v", err)
	}
	if err := os.WriteFile(outputPath, buf.Bytes(), 0o644); err != nil {
		fail("cannot write %s: %v", outputPath, err)
	}
	fmt.Println("Symbolic identities and independent exact matrices passed: r=3..40, k=0,1,2,10,10^6.")
	fmt.Println("Affinely supported real-root witnesses are identical in all checked ranks.")
	fmt.Println("Source SHA256:", data.SourceSHA256)
}
