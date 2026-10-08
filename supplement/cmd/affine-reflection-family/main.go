// Command affine-reflection-family certifies an exact infinite family of
// full-support Gern-terminal elements of the affine Coxeter group of type
// E8 (nine generators, Cartan matrix of the extended E8 diagram).
//
// For beta = (1,2,3,3,2,2,1,1,2) and the null root delta = (2,4,6,5,4,3,2,1,3)
// it certifies, by exact integer arithmetic on the geometric representation:
//   - the Cartan pairing of beta with the simple roots, <beta,beta> = 2 and
//     <delta,e_j> = 0;
//   - an explicit conjugating translation T with T(beta) = beta + delta and
//     T(delta) = delta, obtained from a reflection word for e_3 + delta;
//   - the symbolic inversion count of r_{beta+k delta}: slope 58 and offset
//     -25 from the 240 finite roots (only the finite roots are enumerated,
//     never the affine group);
//   - for k = 0..10, that the reflection r_{beta+k delta} has length
//     33 + 58k, is terminal (left and right), has a reduced word that
//     reproduces the element, and that both one-generator right extensions
//     (generators 0 and 1) are terminal of length 34 + 58k with a
//     full-support fully commutative bottom element below them.
//
// The matrices are 9x9 int64 arrays with an overflow guard on every
// addition and multiplication.
//
// Ports: research/en_families/affine_reflection_family.py, together with the
// part of research/en_families/targeted.py it uses (class Coxeter for en(9):
// right action, descents, reduced word, length, inverse, fully commutative
// test, Bruhat comparison, terminal test).
// Imports: standard library only; no other engine and no internal package.
//
// Usage (inside the work directory): affine-reflection-family
// Writes research/en_families/affine_reflection_family.json and prints the
// same JSON on standard output.
package main

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
)

const n = 9

type vec [n]int64
type mat [n]vec // row i is the image data of generator i, as in the original

var delta = vec{2, 4, 6, 5, 4, 3, 2, 1, 3}
var beta = vec{1, 2, 3, 3, 2, 2, 1, 1, 2}

const outPath = "research/en_families/affine_reflection_family.json"

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "affine-reflection-family: "+format+"\n", args...)
	os.Exit(1)
}

func check(cond bool, what string) {
	if !cond {
		fail("assertion failed: %s", what)
	}
}

const guard = int64(1) << 50

func chk(x int64) int64 {
	if x > guard || x < -guard {
		panic("integer overflow guard tripped")
	}
	return x
}

// coxeter is the simply laced Coxeter group of the diagram en(9): the path
// 0-1-...-7 with node 8 attached to node 2.
type coxeter struct {
	adj      [n][]int
	e        mat
	lengthC  map[mat]int
	descC    map[mat][]int
	fcC      map[mat]bool
	isNeigh  [n][n]bool
	rightMem map[rightKey]mat
}

type rightKey struct {
	w mat
	s int
}

func newCoxeter() *coxeter {
	g := &coxeter{lengthC: map[mat]int{}, descC: map[mat][]int{}, fcC: map[mat]bool{}, rightMem: map[rightKey]mat{}}
	var edges [][2]int
	for i := 0; i < n-2; i++ {
		edges = append(edges, [2]int{i, i + 1})
	}
	edges = append(edges, [2]int{2, n - 1})
	for _, ed := range edges {
		s, t := ed[0], ed[1]
		g.adj[s] = append(g.adj[s], t)
		g.adj[t] = append(g.adj[t], s)
		g.isNeigh[s][t] = true
		g.isNeigh[t][s] = true
	}
	for j := 0; j < n; j++ {
		g.e[j][j] = 1
	}
	return g
}

// clearCaches mirrors the cache_clear calls of the original (memory only;
// results do not depend on it).
func (g *coxeter) clearCaches() {
	g.rightMem = map[rightKey]mat{}
	g.descC = map[mat][]int{}
	g.lengthC = map[mat]int{}
}

func (g *coxeter) right(w mat, s int) mat {
	k := rightKey{w, s}
	if r, ok := g.rightMem[k]; ok {
		return r
	}
	out := w
	for i := 0; i < n; i++ {
		out[s][i] = -w[s][i]
	}
	for _, t := range g.adj[s] {
		for i := 0; i < n; i++ {
			out[t][i] = chk(w[s][i] + w[t][i])
		}
	}
	g.rightMem[k] = out
	return out
}

func (g *coxeter) desc(w mat) []int {
	if d, ok := g.descC[w]; ok {
		return d
	}
	d := []int{}
	for i := 0; i < n; i++ {
		all := true
		for _, x := range w[i] {
			if x > 0 {
				all = false
				break
			}
		}
		if all {
			d = append(d, i)
		}
	}
	g.descC[w] = d
	return d
}

func contains(l []int, x int) bool {
	for _, v := range l {
		if v == x {
			return true
		}
	}
	return false
}

// word returns the reduced word obtained by peeling the smallest right
// descent repeatedly: word(w) = word(w s) + (s,).
func (g *coxeter) word(w mat) []int {
	var rev []int
	for w != g.e {
		ds := g.desc(w)
		check(len(ds) > 0, "non-identity element has a descent")
		s := ds[0]
		rev = append(rev, s)
		w = g.right(w, s)
	}
	out := make([]int, len(rev))
	for i, s := range rev {
		out[len(rev)-1-i] = s
	}
	return out
}

func (g *coxeter) length(w mat) int {
	if l, ok := g.lengthC[w]; ok {
		return l
	}
	l := 0
	cur := w
	var path []mat
	for cur != g.e {
		if v, ok := g.lengthC[cur]; ok {
			l = v
			break
		}
		path = append(path, cur)
		ds := g.desc(cur)
		check(len(ds) > 0, "non-identity element has a descent")
		cur = g.right(cur, ds[0])
	}
	for i := len(path) - 1; i >= 0; i-- {
		l++
		g.lengthC[path[i]] = l
	}
	if w == g.e {
		return 0
	}
	return g.lengthC[w]
}

// elt multiplies the generators of word onto the identity, in order.
func (g *coxeter) elt(word []int) mat {
	w := g.e
	for _, s := range word {
		w = g.right(w, s)
	}
	return w
}

func (g *coxeter) inv(w mat) mat {
	word := g.word(w)
	r := g.e
	for i := len(word) - 1; i >= 0; i-- {
		r = g.right(r, word[i])
	}
	return r
}

// fc reports whether w is fully commutative (no braid of length 3 occurs in
// any reduced word), by the descent criterion of the original.
func (g *coxeter) fc(w mat) bool {
	if v, ok := g.fcC[w]; ok {
		return v
	}
	ds := g.desc(w)
	res := true
	for _, s := range ds {
		for _, t := range g.adj[s] {
			if contains(ds, t) {
				res = false
			}
		}
	}
	if res {
		for _, s := range ds {
			if !g.fc(g.right(w, s)) {
				res = false
				break
			}
		}
	}
	g.fcC[w] = res
	return res
}

// terminal reports the Gern-terminal property: for v = w and v = w^-1, no
// right descent s of v has a neighbour t that becomes a descent of v s.
func (g *coxeter) terminal(w mat) bool {
	for _, v := range []mat{w, g.inv(w)} {
		for _, s := range g.desc(v) {
			vs := g.desc(g.right(v, s))
			for _, t := range g.adj[s] {
				if contains(vs, t) {
					return false
				}
			}
		}
	}
	return true
}

// leq is the Bruhat comparison x <= w by the descent recursion of the
// original.
func (g *coxeter) leq(x, w mat) bool {
	for x != w {
		if g.length(x) >= g.length(w) {
			return false
		}
		s := g.desc(w)[0]
		if contains(g.desc(x), s) {
			x = g.right(x, s)
		}
		w = g.right(w, s)
	}
	return true
}

// pairing is the symmetric bilinear form of the Cartan matrix.
func (g *coxeter) pairing(a, b vec) int64 {
	var total int64
	for i := 0; i < n; i++ {
		inner := 2 * b[i]
		for _, j := range g.adj[i] {
			inner -= b[j]
		}
		total = chk(total + chk(a[i]*inner))
	}
	return total
}

func addScaled(a, b vec, k int64) vec {
	var c vec
	for i := range a {
		c[i] = chk(a[i] + chk(k*b[i]))
	}
	return c
}

// action is the matrix action: (w a)_i = sum_j w[j][i] a[j].
func action(w mat, a vec) vec {
	var out vec
	for i := 0; i < n; i++ {
		var t int64
		for j := 0; j < n; j++ {
			t = chk(t + chk(w[j][i]*a[j]))
		}
		out[i] = t
	}
	return out
}

// rootReflection is the reflection r_a as a matrix (row j = e_j - <e_j,a> a).
func (g *coxeter) rootReflection(a vec) mat {
	var out mat
	for j := 0; j < n; j++ {
		out[j] = addScaled(g.e[j], a, -g.pairing(g.e[j], a))
	}
	return out
}

func (g *coxeter) reflect(a vec, s int) vec {
	out := a
	out[s] = chk(out[s] - g.pairing(g.e[s], a))
	return out
}

func minOf(a vec) int64 {
	m := a[0]
	for _, v := range a {
		if v < m {
			m = v
		}
	}
	return m
}

func sumOf(a vec) int64 {
	var t int64
	for _, v := range a {
		t += v
	}
	return t
}

// reflectionWord returns a palindromic word for the reflection of the
// positive real root a, peeling positive-pairing simple roots.
func (g *coxeter) reflectionWord(a vec) []int {
	for i := 0; i < n; i++ {
		if a == g.e[i] {
			return []int{i}
		}
	}
	s := -1
	for i := 0; i < n; i++ {
		if g.pairing(g.e[i], a) > 0 {
			s = i
			break
		}
	}
	check(s >= 0, "root has a simple root of positive pairing")
	b := g.reflect(a, s)
	check(minOf(b) >= 0 && sumOf(b) < sumOf(a), "reflection step lowers the root")
	inner := g.reflectionWord(b)
	out := make([]int, 0, len(inner)+2)
	out = append(out, s)
	out = append(out, inner...)
	out = append(out, s)
	return out
}

// finiteRoots enumerates the 240 roots of the parabolic subgroup obtained by
// deleting node 7, sorted lexicographically.
func (g *coxeter) finiteRoots() []vec {
	var inds []int
	for i := 0; i < n; i++ {
		if i != 7 {
			inds = append(inds, i)
		}
	}
	seen := map[vec]bool{}
	var todo []vec
	for _, i := range inds {
		seen[g.e[i]] = true
		todo = append(todo, g.e[i])
	}
	for idx := 0; idx < len(todo); idx++ {
		a := todo[idx]
		for _, s := range inds {
			b := g.reflect(a, s)
			if !seen[b] {
				seen[b] = true
				todo = append(todo, b)
			}
		}
	}
	roots := append([]vec(nil), todo...)
	sort.Slice(roots, func(i, j int) bool {
		for k := 0; k < n; k++ {
			if roots[i][k] != roots[j][k] {
				return roots[i][k] < roots[j][k]
			}
		}
		return false
	})
	return roots
}

type rootClass struct {
	Pairing          int64 `json:"pairing"`
	MinimumAffineLvl int   `json:"minimum_affine_level"`
	ImageNegative    int   `json:"image_negative"`
	Count            int   `json:"count"`
}

type extension struct {
	Generator    int   `json:"generator"`
	Length       int   `json:"length"`
	Desc         []int `json:"desc"`
	BottomWord   []int `json:"bottom_word"`
	BottomLength int   `json:"bottom_length"`
	IntervalRank int   `json:"interval_rank"`
}

type row struct {
	K          int         `json:"k"`
	Beta       []int64     `json:"beta"`
	Length     int         `json:"length"`
	Desc       []int       `json:"desc"`
	Terminal   bool        `json:"terminal"`
	Extensions []extension `json:"extensions"`
}

type output struct {
	Delta                  []int64     `json:"delta"`
	Beta                   []int64     `json:"beta"`
	CartanPairing          []int64     `json:"cartan_pairing"`
	BaseReflectionWord     []int       `json:"base_reflection_word"`
	TranslationWord        []int       `json:"translation_word"`
	TranslationReducedWord []int       `json:"translation_reduced_word"`
	FiniteRootCount        int         `json:"finite_root_count"`
	LengthSlope            int64       `json:"length_slope"`
	LengthOffsetForKPlus1  int64       `json:"length_offset_for_k_plus_1"`
	InversionRootClasses   []rootClass `json:"inversion_root_classes"`
	Checks                 []row       `json:"checks"`
}

func slice(a vec) []int64 { return append([]int64(nil), a[:]...) }

func build() output {
	g := newCoxeter()

	var m vec
	for j := 0; j < n; j++ {
		m[j] = g.pairing(g.e[j], beta)
	}
	check(m == vec{0, 0, -1, 1, -1, 1, -1, 1, 1}, "Cartan pairing of beta with the simple roots")
	for j := 0; j < n; j++ {
		check(g.pairing(delta, g.e[j]) == 0, "delta is orthogonal to every simple root")
	}
	check(g.pairing(beta, beta) == 2, "beta is a real root")

	// Root membership and an explicit conjugating translation certificate.
	a := addScaled(g.e[3], delta, 1)
	rword := g.reflectionWord(a)
	check(g.elt(rword) == g.rootReflection(a), "reflection word of e3+delta")
	translationWord := append(append([]int(nil), rword...), 3)
	T := g.elt(translationWord)
	check(action(T, beta) == addScaled(beta, delta, 1), "T(beta) = beta + delta")
	check(action(T, delta) == delta, "T(delta) = delta")
	bword := g.reflectionWord(beta)
	check(g.elt(bword) == g.rootReflection(beta), "reflection word of beta")

	// Finite-root counting gives the affine inversion count symbolically.
	b0 := addScaled(beta, delta, -1)
	roots := g.finiteRoots()
	check(len(roots) == 240, "240 finite roots")
	rootSet := map[vec]bool{}
	for _, r := range roots {
		rootSet[r] = true
	}
	var slope, offset int64
	type ckey struct {
		val        int64
		nmin, imgn int
	}
	counts := map[ckey]int{}
	for _, r := range roots {
		val := g.pairing(r, b0)
		if val <= 0 {
			continue
		}
		image := addScaled(r, b0, -val)
		check(rootSet[image], "image of a root lies in the finite root system")
		nmin := 0
		if minOf(r) < 0 {
			nmin = 1
		}
		imageNeg := 0
		if minOf(image) < 0 {
			imageNeg = 1
		}
		slope += val
		offset += int64(-nmin + imageNeg)
		counts[ckey{val, nmin, imageNeg}]++
	}
	check(slope == 58 && offset == -25, "slope 58 and offset -25")
	keys := make([]ckey, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		a, b := keys[i], keys[j]
		if a.val != b.val {
			return a.val < b.val
		}
		if a.nmin != b.nmin {
			return a.nmin < b.nmin
		}
		return a.imgn < b.imgn
	})
	classes := []rootClass{}
	for _, k := range keys {
		classes = append(classes, rootClass{k.val, k.nmin, k.imgn, counts[k]})
	}

	checks := []row{}
	for k := 0; k <= 10; k++ {
		bk := addScaled(beta, delta, int64(k))
		w := g.rootReflection(bk)
		check(g.length(w) == 33+58*k && g.terminal(w), fmt.Sprintf("k=%d: length 33+58k and terminal", k))
		check(g.elt(g.word(w)) == w, fmt.Sprintf("k=%d: reduced word reproduces the element", k))
		r := row{K: k, Beta: slice(bk), Length: g.length(w), Desc: g.desc(w), Terminal: g.terminal(w), Extensions: []extension{}}
		for s := 0; s <= 1; s++ {
			v := g.right(w, s)
			check(g.length(v) == 34+58*k && g.terminal(v), fmt.Sprintf("k=%d s=%d: extension has length 34+58k and is terminal", k, s))
			x := g.elt(g.desc(v))
			check(g.fc(x) && g.leq(x, v), fmt.Sprintf("k=%d s=%d: bottom element is fully commutative and below the extension", k, s))
			r.Extensions = append(r.Extensions, extension{
				Generator: s, Length: g.length(v), Desc: g.desc(v), BottomWord: g.desc(v),
				BottomLength: g.length(x), IntervalRank: g.length(v) - g.length(x)})
		}
		checks = append(checks, r)
		g.clearCaches()
	}

	baseWord := g.word(g.rootReflection(beta))
	reduced := g.word(T)
	return output{
		Delta: slice(delta), Beta: slice(beta), CartanPairing: slice(m),
		BaseReflectionWord: baseWord, TranslationWord: translationWord,
		TranslationReducedWord: reduced, FiniteRootCount: len(roots),
		LengthSlope: slope, LengthOffsetForKPlus1: offset,
		InversionRootClasses: classes, Checks: checks,
	}
}

func main() {
	out := build()
	data, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		fail("%v", err)
	}
	data = append(data, '\n')
	if err := os.MkdirAll(filepath.Dir(outPath), 0o755); err != nil {
		fail("%v", err)
	}
	if err := os.WriteFile(outPath, data, 0o644); err != nil {
		fail("%v", err)
	}
	os.Stdout.Write(data)
}
