// Command cartan-candidates exhausts the terminal reflections r_beta of the
// Coxeter group W(E_n) (n >= 10 in the committed runs) that have a prescribed
// set of small root pairings with the simple roots, and certifies each one it
// finds.
//
// For a reflection r_beta (beta a norm-2 vector, i.e. a real root) the
// pairings m_j = (beta, alpha_j) determine beta = A^{-1} m, where A is the
// Cartan matrix.  The search enumerates
//   - a nonempty independent set I of the Dynkin diagram (all of them, or only
//     those of maximal size with -max-only),
//   - the pairings m_s in 1..M on I (these are the right and left descents of
//     r_beta), and
//   - the pairings m_j in -M..-demand(j) on the other nodes, where demand(j)
//     is the largest m_s on the neighbours s of j that lie in I,
//
// keeps the vectors beta = A^{-1} m that are positive integer vectors of norm
// 2 ((beta, m) = 2), certifies beta as a real root by reducing it to a simple
// root with simple reflections (height strictly decreasing), builds the
// reflection as an integer matrix and asserts that it is terminal, an
// involution and has descent set exactly I, and that the word conjugator +
// simple root + reversed conjugator reproduces it.  All arithmetic is exact
// (big.Rat for the inverse Cartan matrix; int64 with overflow guards for the
// rest).
//
// Ports research/en_families/cartan_candidates.py (with the Coxeter-group
// arithmetic of research/en_families/targeted.py that it uses: right
// multiplication on integer matrices, descent sets, the canonical reduced
// word, inverses, the terminal test).  Standard library only; it imports no
// other package of this module (self-contained).  E_n labelling: chain
// 0-1-...-(n-2), node n-1 attached to node 2.  Elements are integer matrices
// whose row i is the image data of node i exactly as in the original.
//
// Usage: cartan-candidates -rank N [-max-entry M] [-max-only] [-out FILE]
// Writes research/en_families/cartan_E{N}_m{M}_max{0|1}.json (relative to the
// current directory) unless -out is given; the original's "seconds" key is
// dropped, everything else is unchanged.  Each root found is also printed as a
// JSON line on stdout, followed by a SUMMARY line.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/big"
	"os"
)

// ---------- checked int64 arithmetic ----------

func addC(a, b int64) int64 {
	c := a + b
	if (c > a) != (b > 0) {
		panic("int64 overflow in addition")
	}
	return c
}

func subC(a, b int64) int64 {
	c := a - b
	if (c < a) != (b > 0) {
		panic("int64 overflow in subtraction")
	}
	return c
}

func mulC(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	c := a * b
	if c/b != a || (a == math.MinInt64 && b == -1) || (b == math.MinInt64 && a == -1) {
		panic("int64 overflow in multiplication")
	}
	return c
}

// ---------- the Coxeter group E_n (simply laced, integer matrices) ----------

type group struct {
	n   int
	adj [][]int
	e   []int64 // identity matrix, row-major
}

type elt []int64 // n*n row-major; row i is e[i*n:(i+1)*n]

func en(n int) *group {
	if n < 4 {
		fatalf("rank must be at least 4")
	}
	g := &group{n: n, adj: make([][]int, n)}
	edge := func(a, b int) {
		g.adj[a] = append(g.adj[a], b)
		g.adj[b] = append(g.adj[b], a)
	}
	for i := 0; i < n-2; i++ {
		edge(i, i+1)
	}
	edge(2, n-1)
	g.e = make([]int64, n*n)
	for i := 0; i < n; i++ {
		g.e[i*n+i] = 1
	}
	return g
}

func (g *group) isAdj(s, t int) bool {
	for _, x := range g.adj[s] {
		if x == t {
			return true
		}
	}
	return false
}

// right returns w*s: row s is negated and row t becomes w[s]+w[t] for every
// neighbour t of s.
func (g *group) right(w elt, s int) elt {
	n := g.n
	out := make(elt, len(w))
	copy(out, w)
	for k := 0; k < n; k++ {
		out[s*n+k] = -w[s*n+k]
		if w[s*n+k] == math.MinInt64 {
			panic("int64 overflow in negation")
		}
	}
	for _, t := range g.adj[s] {
		for k := 0; k < n; k++ {
			out[t*n+k] = addC(w[s*n+k], w[t*n+k])
		}
	}
	return out
}

// desc returns the ascending list of rows that are entirely <= 0.
func (g *group) desc(w elt) []int {
	n := g.n
	var out []int
	for i := 0; i < n; i++ {
		ok := true
		for k := 0; k < n; k++ {
			if w[i*n+k] > 0 {
				ok = false
				break
			}
		}
		if ok {
			out = append(out, i)
		}
	}
	return out
}

func (g *group) isIdentity(w elt) bool {
	for i, v := range w {
		if v != g.e[i] {
			return false
		}
	}
	return true
}

// word returns the canonical reduced word of w: word(w) = word(w s) + (s,)
// with s the smallest right descent.
func (g *group) word(w elt) []int {
	var rev []int
	for !g.isIdentity(w) {
		ds := g.desc(w)
		if len(ds) == 0 {
			fatalf("element has no right descent but is not the identity")
		}
		rev = append(rev, ds[0])
		w = g.right(w, ds[0])
	}
	out := make([]int, len(rev))
	for i, s := range rev {
		out[len(rev)-1-i] = s
	}
	return out
}

func (g *group) elt(word []int) elt {
	w := make(elt, len(g.e))
	copy(w, g.e)
	for _, s := range word {
		w = g.right(w, s)
	}
	return w
}

func (g *group) inv(w elt) elt {
	word := g.word(w)
	rev := make([]int, len(word))
	for i, s := range word {
		rev[len(word)-1-i] = s
	}
	return g.elt(rev)
}

func equalElt(a, b elt) bool {
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

// terminal: for v in (w, w^-1) and every descent s of v, no neighbour of s
// is a descent of v s.
func (g *group) terminal(w elt) bool {
	for _, v := range []elt{w, g.inv(w)} {
		for _, s := range g.desc(v) {
			ds := g.desc(g.right(v, s))
			for _, t := range g.adj[s] {
				for _, d := range ds {
					if d == t {
						return false
					}
				}
			}
		}
	}
	return true
}

// ---------- the search ----------

// inverseCartan returns the inverse of the Cartan matrix as exact rationals
// (Gauss-Jordan elimination, first nonzero pivot).
func inverseCartan(g *group) [][]*big.Rat {
	n := g.n
	A := make([][]*big.Rat, n)
	for i := 0; i < n; i++ {
		A[i] = make([]*big.Rat, 2*n)
		for j := 0; j < n; j++ {
			v := int64(0)
			if i == j {
				v = 2
			} else if g.isAdj(i, j) {
				v = -1
			}
			A[i][j] = new(big.Rat).SetInt64(v)
			id := int64(0)
			if i == j {
				id = 1
			}
			A[i][n+j] = new(big.Rat).SetInt64(id)
		}
	}
	for j := 0; j < n; j++ {
		pivot := -1
		for i := j; i < n; i++ {
			if A[i][j].Sign() != 0 {
				pivot = i
				break
			}
		}
		if pivot < 0 {
			fatalf("singular Cartan matrix")
		}
		A[j], A[pivot] = A[pivot], A[j]
		v := new(big.Rat).Set(A[j][j])
		for k := range A[j] {
			A[j][k] = new(big.Rat).Quo(A[j][k], v)
		}
		for i := 0; i < n; i++ {
			if i == j {
				continue
			}
			v := new(big.Rat).Set(A[i][j])
			if v.Sign() != 0 {
				for k := range A[i] {
					A[i][k] = new(big.Rat).Sub(A[i][k], new(big.Rat).Mul(v, A[j][k]))
				}
			}
		}
	}
	out := make([][]*big.Rat, n)
	for i := range out {
		out[i] = A[i][n:]
	}
	return out
}

// independentSets lists the independent sets of the Dynkin diagram in the
// order of the original recursion (exclude node i first, then include it).
func independentSets(g *group) [][]int {
	var out [][]int
	var visit func(i int, chosen []int)
	visit = func(i int, chosen []int) {
		if i == g.n {
			out = append(out, append([]int{}, chosen...))
			return
		}
		visit(i+1, chosen)
		ok := true
		for _, j := range g.adj[i] {
			for _, c := range chosen {
				if c == j {
					ok = false
				}
			}
		}
		if ok {
			visit(i+1, append(append([]int{}, chosen...), i))
		}
	}
	visit(0, nil)
	return out
}

type rootReduction struct {
	Conjugator []int `json:"conjugator"`
	SimpleRoot int   `json:"simple_root"`
}

// reduceRoot reduces the positive vector b to a simple root by simple
// reflections that strictly lower the height; nil if it fails.
func reduceRoot(g *group, b []int64) *rootReduction {
	n := g.n
	steps := []int{}
	a := append([]int64{}, b...)
	sum := func(v []int64) int64 {
		var t int64
		for _, x := range v {
			t = addC(t, x)
		}
		return t
	}
	for sum(a) > 1 {
		s := -1
		var ps int64
		for i := 0; i < n; i++ {
			p := mulC(2, a[i])
			for _, j := range g.adj[i] {
				p = subC(p, a[j])
			}
			if p > 0 {
				s, ps = i, p
				break
			}
		}
		if s < 0 {
			return nil
		}
		v := append([]int64{}, a...)
		v[s] = subC(v[s], ps)
		for _, x := range v {
			if x < 0 {
				return nil
			}
		}
		a = v
		steps = append(steps, s)
	}
	for i := 0; i < n; i++ {
		isUnit := true
		for k := 0; k < n; k++ {
			want := int64(0)
			if k == i {
				want = 1
			}
			if a[k] != want {
				isUnit = false
				break
			}
		}
		if isUnit {
			return &rootReduction{Conjugator: steps, SimpleRoot: i}
		}
	}
	return nil
}

type rootRecord struct {
	Rank               int            `json:"rank"`
	Beta               []int64        `json:"beta"`
	Pairings           []int64        `json:"pairings"`
	I                  []int          `json:"I"`
	MaximumIndependent bool           `json:"maximum_independent"`
	IndependenceNumber int            `json:"independence_number"`
	Height             int64          `json:"height"`
	Length             int            `json:"length"`
	RootReduction      *rootReduction `json:"root_reduction"`
	ReducedWord        []int          `json:"reduced_word"`
}

type result struct {
	Rank                   int          `json:"rank"`
	MaxPairingAbsoluteVal  int          `json:"max_pairing_absolute_value"`
	MaximumIndependentOnly bool         `json:"maximum_independent_only"`
	IndependenceNumber     int          `json:"independence_number"`
	IndependentSetsTested  int          `json:"independent_sets_tested"`
	PairingsTested         int          `json:"pairings_tested"`
	PositiveIntegerVectors int          `json:"positive_integer_vectors"`
	NormTwoVectors         int          `json:"norm_two_vectors"`
	RealTerminalRoots      []rootRecord `json:"real_terminal_roots"`
}

func gcd(a, b *big.Int) *big.Int { return new(big.Int).GCD(nil, nil, a, b) }

func search(n, maxEntry int, maxOnly bool, emit func(rootRecord)) result {
	g := en(n)
	inv := inverseCartan(g)
	denom := big.NewInt(1)
	for _, row := range inv {
		for _, v := range row {
			d := v.Denom()
			l := new(big.Int).Mul(denom, d)
			denom = l.Quo(l, gcd(denom, d))
		}
	}
	if !denom.IsInt64() {
		panic("common denominator exceeds int64")
	}
	den := denom.Int64()
	N := make([][]int64, n)
	for i, row := range inv {
		N[i] = make([]int64, n)
		for j, v := range row {
			x := new(big.Rat).Mul(v, new(big.Rat).SetInt(denom))
			if !x.IsInt() || !x.Num().IsInt64() {
				panic("scaled inverse Cartan entry is not an int64")
			}
			N[i][j] = x.Num().Int64()
		}
	}
	all := independentSets(g)
	alpha := 0
	for _, I := range all {
		if len(I) > alpha {
			alpha = len(I)
		}
	}
	var Is [][]int
	for _, I := range all {
		if len(I) > 0 && (!maxOnly || len(I) == alpha) {
			Is = append(Is, I)
		}
	}
	res := result{Rank: n, MaxPairingAbsoluteVal: maxEntry, MaximumIndependentOnly: maxOnly,
		IndependenceNumber: alpha, IndependentSetsTested: len(Is), RealTerminalRoots: []rootRecord{}}
	inI := make([]bool, n)
	m := make([]int64, n)
	for _, I := range Is {
		for k := range inI {
			inI[k] = false
		}
		for _, s := range I {
			inI[s] = true
		}
		var rest []int
		for j := 0; j < n; j++ {
			if !inI[j] {
				rest = append(rest, j)
			}
		}
		pos := make([]int, len(I)) // 0-based digits of the pairings on I (value digit+1)
		for {
			for k, s := range I {
				m[s] = int64(pos[k] + 1)
			}
			// ranges for the other nodes
			lo := make([]int64, len(rest))
			hi := make([]int64, len(rest))
			for k, j := range rest {
				var demand int64
				for _, s := range g.adj[j] {
					if inI[s] && m[s] > demand {
						demand = m[s]
					}
				}
				lo[k] = -int64(maxEntry)
				hi[k] = -demand
			}
			vals := append([]int64{}, lo...)
			for {
				res.PairingsTested++
				for k, j := range rest {
					m[j] = vals[k]
				}
				processCandidate(g, &res, I, m, N, den, alpha, emit)
				// odometer: last position varies fastest
				k := len(rest) - 1
				for ; k >= 0; k-- {
					if vals[k] < hi[k] {
						vals[k]++
						break
					}
					vals[k] = lo[k]
				}
				if k < 0 {
					break
				}
			}
			k := len(I) - 1
			for ; k >= 0; k-- {
				if pos[k] < maxEntry-1 {
					pos[k]++
					break
				}
				pos[k] = 0
			}
			if k < 0 {
				break
			}
		}
	}
	return res
}

func processCandidate(g *group, res *result, I []int, m []int64, N [][]int64, den int64, alpha int, emit func(rootRecord)) {
	n := g.n
	b := make([]int64, n)
	for i := 0; i < n; i++ {
		var x int64
		for k := 0; k < n; k++ {
			x = addC(x, mulC(N[i][k], m[k]))
		}
		if x <= 0 || x%den != 0 {
			return
		}
		b[i] = x / den
	}
	res.PositiveIntegerVectors++
	var norm int64
	for i := 0; i < n; i++ {
		norm = addC(norm, mulC(b[i], m[i]))
	}
	if norm != 2 {
		return
	}
	res.NormTwoVectors++
	cert := reduceRoot(g, b)
	if cert == nil {
		return
	}
	// w = rows e_j - m_j * b
	w := make(elt, n*n)
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			e := int64(0)
			if i == j {
				e = 1
			}
			w[j*n+i] = subC(e, mulC(m[j], b[i]))
		}
	}
	if !g.terminal(w) {
		fatalf("assertion failed: reflection for beta=%v is not terminal", b)
	}
	if !equalInts(g.desc(w), I) {
		fatalf("assertion failed: descent set of reflection for beta=%v differs from I=%v", b, I)
	}
	if !equalElt(g.inv(w), w) {
		fatalf("assertion failed: reflection for beta=%v is not an involution", b)
	}
	c := cert.Conjugator
	word := append([]int{}, c...)
	word = append(word, cert.SimpleRoot)
	for i := len(c) - 1; i >= 0; i-- {
		word = append(word, c[i])
	}
	if !equalElt(g.elt(word), w) {
		fatalf("assertion failed: conjugator word does not reproduce the reflection for beta=%v", b)
	}
	rw := g.word(w)
	var height int64
	for _, x := range b {
		height = addC(height, x)
	}
	rec := rootRecord{Rank: n, Beta: b, Pairings: append([]int64{}, m...), I: append([]int{}, I...),
		MaximumIndependent: len(I) == alpha, IndependenceNumber: alpha, Height: height,
		Length: len(rw), RootReduction: cert, ReducedWord: rw}
	res.RealTerminalRoots = append(res.RealTerminalRoots, rec)
	if emit != nil {
		emit(rec)
	}
}

func fatalf(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "cartan-candidates: "+format+"\n", args...)
	os.Exit(1)
}

func main() {
	rank := flag.Int("rank", 0, "rank n of E_n (required, >= 4)")
	maxEntry := flag.Int("max-entry", 1, "bound M on the absolute value of the root pairings")
	maxOnly := flag.Bool("max-only", false, "only maximum-size independent sets I")
	out := flag.String("out", "", "output file (default research/en_families/cartan_E{N}_m{M}_max{0|1}.json)")
	flag.Parse()
	if flag.NArg() != 0 || *rank < 4 || *maxEntry < 1 {
		fmt.Fprintln(os.Stderr, "usage: cartan-candidates -rank N [-max-entry M] [-max-only] [-out FILE]")
		os.Exit(2)
	}
	res := search(*rank, *maxEntry, *maxOnly, func(r rootRecord) {
		b, err := json.Marshal(r)
		if err != nil {
			fatalf("%v", err)
		}
		fmt.Println(string(b))
	})
	path := *out
	if path == "" {
		mo := 0
		if *maxOnly {
			mo = 1
		}
		path = fmt.Sprintf("research/en_families/cartan_E%d_m%d_max%d.json", *rank, *maxEntry, mo)
	}
	data, err := json.MarshalIndent(res, "", "  ")
	if err != nil {
		fatalf("%v", err)
	}
	data = append(data, '\n')
	if err := os.WriteFile(path, data, 0o644); err != nil {
		fatalf("%v", err)
	}
	fmt.Printf("SUMMARY rank=%d max_pairing_absolute_value=%d maximum_independent_only=%t independence_number=%d independent_sets_tested=%d pairings_tested=%d positive_integer_vectors=%d norm_two_vectors=%d roots=%d\n",
		res.Rank, res.MaxPairingAbsoluteVal, res.MaximumIndependentOnly, res.IndependenceNumber,
		res.IndependentSetsTested, res.PairingsTested, res.PositiveIntegerVectors, res.NormTwoVectors,
		len(res.RealTerminalRoots))
}
